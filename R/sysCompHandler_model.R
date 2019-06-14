#  nucdataBaynet - Nuclear Data Evaluation Using a Bayesian Network 
#  Copyright (C) 2019  Georg Schnabel
#  
#  nucdataBaynet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  nucdataBaynet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>


#' Create Handler for Systematic Model Components
#'
#' @return TODO
#' @import data.table Matrix methods
#' @export
createSysCompModelHandler <- function() {

  refParamDt <- NULL
  refNeedsDt <- NULL
  refSparDt <- NULL
  priorParamDt <- NULL
  exforHandler <- NULL
  subentList <- NULL
  thisTrafo <- NULL
  # for caching
  saveSexp <- NULL
  saveExpDt <- NULL

  # administration

  setRef <- function(needsDt, SparDt, paramDt, exforHandler, subents) {

    stopifnot(all(c("IDX", "V1") %in% colnames(needsDt)))
    stopifnot(all(c("IDX1", "IDX2", "X") %in% colnames(SparDt)))
    stopifnot(all(c("PROJECTILE", "ELEMENT", "MASS", "PARNAME", "PARVAL",
                    "ADJUSTABLE", "IDX") %in% colnames(paramDt)))

    adjParamDt <- copy(paramDt)
    setkey(adjParamDt, IDX)
    adjParamDt <- adjParamDt[ADJUSTABLE == TRUE]
    adjParamDt[, NEWIDX := seq_len(.N)]
    setkey(adjParamDt, IDX)
    adjSparDt <- copy(SparDt)
    adjSparDt[, IDX2 := adjParamDt[J(IDX2), NEWIDX]]
    stopifnot(all(!is.na(adjSparDt$IDX2)))
    adjParamDt[, IDX := NEWIDX]
    adjParamDt[, NEWIDX := NULL]
    adjParamDt[, PARVAL := unlist(PARVAL)]
    stopifnot(is.numeric(adjParamDt$PARVAL))

    refNeedsDt <<- needsDt
    refParamDt <<- adjParamDt
    refSparDt <<- adjSparDt
    exforHandler <<- exforHandler
    subentList <<- subents
    saveSexp <<- NULL
    saveExpDt <<- NULL
  }

  setPrior <- function(paramDt) {

    stopifnot(all(c("PROJECTILE", "ELEMENT", "MASS", "PARNAME", "PARVAL",
                    "PARUNC", "ADJUSTABLE", "IDX") %in% colnames(paramDt)))
    newPriorParamDt <<- copy(refParamDt)
    newPriorParamDt[, PARVAL := NA_real_]
    setkey(newPriorParamDt, PARNAME)
    newPriorParamDt[J(paramDt[ADJUSTABLE == TRUE, PARNAME]),
                 PARVAL := paramDt[ADJUSTABLE == TRUE, unlist(PARVAL)]]
    if (!is.numeric(newPriorParamDt$PARVAL) || any(is.na(newPriorParamDt$PARVAL))) {
      stop(paste0("paramDt does not contain all required parameters, such as ",
                  newPriorParamDt[is.na(PARVAL),
                                  paste0(PARNAME[seq_len(min(5,length(PARNAME)))], collapse = ", ")]))
    }
    priorParamDt <<- newPriorParamDt
    TRUE
  }

  setTrafo <- function(trafo) {

    stopifnot(is.list(trafo))
    stopifnot(is.function(trafo$fun))
    stopifnot(is.function(trafo$jac))
    stopifnot(is.function(trafo$invfun))
    thisTrafo <<- trafo
  }

  getTrafo <- function() {

    thisTrafo
  }

  # must-have functions

  createSysDt <- function() {

    setkey(priorParamDt, IDX)
    setkey(refParamDt, IDX)
    stopifnot(all(priorParamDt$IDX == refParamDt$IDX))
    isGP <- grepl("adjust\\(", priorParamDt$PARNAME)
    parnames <- sub("adjust\\([^)]*\\)", "adjust", priorParamDt$PARNAME)
    enVals <- rep(NA_real_, nrow(priorParamDt))
    enVals[isGP] <- as.numeric(sub("^.*adjust\\(([^)]*)\\).*$", "\\1", priorParamDt$PARNAME[isGP]))

    trafoSuffix <- if (!is.null(thisTrafo)) "_trafo" else ""
    priorParamDt[, list(EXPID = paste0("TALYS-",parnames),
                       DIDX = IDX,
                       ERRTYPE = paste0(ifelse(isGP, "talyspar_endep", "talyspar"),
                                        trafoSuffix),
                       GPTYPE = NA_character_,
                       DATA = as.numeric(unlist(PARVAL)),
                       UNC = PARUNC,
                       REFDATA = refParamDt[, PARVAL],
                       EN = enVals)]
  }

  getSignature <- function() {
    "TALYSmodel"
  }


  getErrTypes <- function() {
    c("talyspar", "talyspar_endep",
      "talyspar_trafo", "talyspar_endep_trafo")
  }

  # mapping

  map <- function(expDt, sysDt, ret.mat) {

    curSparDt <- merge(sysDt[ERRTYPE %in% c("talyspar", "talyspar_endep",
                                            "talyspar_trafo", "talyspar_endep_trafo"),
                         list(IDX2 = DIDX, REALIDX = IDX)],
                   refSparDt, by = "IDX2")
    curSparDt[, IDX2 := NULL]
    setnames(curSparDt, "REALIDX", "IDX2")
    setcolorder(curSparDt, c("IDX1","IDX2","X"))
    curSparDt <- curSparDt[X != 0, ]
    resDims <- c(nrow(refNeedsDt), max(nrow(sysDt), max(sysDt$IDX)))
    curSpar <- returnSparseMatrix(curSparDt, resDims, ret.mat = TRUE)

    # some hack to avoid recalculation of Sexp
    curSexp <- NULL
    if (!is.null(saveExpDt) && nrow(expDt) == nrow(saveExpDt))
      if (nrow(merge(expDt[, list(IDX, EXPID, DIDX)], saveExpDt,
                     by = c("IDX", "EXPID", "DIDX"))) == nrow(expDt))
        curSexp <- saveSexp
    if (is.null(curSexp)) {
      curSexp <- exforHandler$getJac(expDt, refNeedsDt, subentList)
      saveSexp <<- curSexp
      saveExpDt <<- expDt[, list(IDX, EXPID, DIDX)]
    }

    if (!is.null(thisTrafo)) {
      trafoSel <- grepl("talyspar(_endep)?(trafo)", sysDt$ERRTYPE)
      trafoIdcs <- sysDt[trafoSel, IDX]
      trafoRefVals <- sysDt[trafoSel, REFDATA]
      Strafo <- thisTrafo$jac(trafoRefVals)
      curSpar[, trafoIdcs] <- curSpar[, trafoIdcs] %*% Strafo
    }

    curSall <- as(curSexp %*% curSpar, "sparseMatrix")
    curSall <- drop0(curSall)
    if (isTRUE(ret.mat))
      curSall
    else {
      resDt <- as.data.table(summary(curSall))
      setnames(resDt, c("i", "j", "x"), c("IDX1", "IDX2", "X"))
    }
  }

  list(setRef = setRef,
       setPrior = setPrior,
       setTrafo = setTrafo,
       getTrafo = getTrafo,
       createSysDt = createSysDt,
       # must have properties
       getSignature = getSignature,
       getErrTypes = getErrTypes,
       map = map)
}
