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


#' Create Handler for Systematic Normalization Components
#'
#' @return TODO
#' @import data.table Matrix
#' @export
#'
createSysCompNormHandler <- function(dataref = NULL) {

  assocDt <- data.table(expid = character(0),
                        cat = character(0), lev = character(0),
                        val = numeric(0), unc = numeric(0),
                        rel = logical(0))

  addSysUnc <- function(category, levels, vals, uncs, rel = TRUE) {

    stopifnot(is.character(category))
    stopifnot(is.character(levels))
    stopifnot(is.numeric(vals))
    stopifnot(is.numeric(uncs))
    newDt <- data.table(expid = paste0(category, "-", levels),
                        cat = category,
                        lev = levels,
                        val = vals,
                        unc = uncs,
                        rel = rel)
    proposedAssocDt <- rbind(assocDt, newDt)
    if (anyDuplicated(proposedAssocDt[, list(cat, lev)]))
      stop("some new cat/lev combination in conflict with existing specifications")
    else
      assocDt <<- proposedAssocDt
  }


  createSysDt <- function() {

    assocDt[, list(EXPID = expid,
                   DIDX = 1L,
                   ERRTYPE = ifelse(rel, "sys-rel", "sys-abs"),
                   GPTYPE = NA_character_,
                   DATA = val,
                   UNC = unc)]
  }


  getErrTypes <- function() {
    c("sys-abs", "sys-rel")
  }


  getSignature <- function() {
    "normHandler"
  }


  map <- function(expDt, sysDt, ret.mat) {

    stopifnot(all(sysDt$ERRTYPE == "sys-abs" | sysDt$ERRTYPE == "sys-rel"))
    stopifnot("IDX" %in% colnames(expDt),
              "IDX" %in% colnames(sysDt))
    if (any(sysDt$ERRTYPE == "sys-rel") && ! dataref %in% names(expDt))
      stop(paste0("Column ", dataref, " specified as reference data ",
                  "must be present in expDt if any ERRTYPE == sys-abs"))
    
    uniqueCats <- unique(assocDt$cat)

    moltenExpDt <- melt(expDt, id.vars = c("IDX", dataref), measure.vars = uniqueCats,
                     variable.name = "CAT", value.name = "LEV")
    moltenExpDt[, EXPID := paste(CAT, LEV, sep="-")]
    moltenExpDt[, DATAREF := get(dataref)]
    setnames(moltenExpDt, "IDX", "IDX1")

    resDt <- merge(moltenExpDt, sysDt, by = "EXPID")
    setnames(resDt, "IDX", "IDX2")
    resDt <- resDt[, list(IDX1, IDX2, ERRTYPE, DATAREF)]
    resDt[, X := ifelse(ERRTYPE == "sys-abs", 1, DATAREF)]
    resDt[, DATAREF := NULL]
    resDt[, ERRTYPE := NULL]

    resDims <- c(nrow(expDt), nrow(sysDt))
    returnSparseMatrix(resDt, resDims, ret.mat)
  }


  list(addSysUnc = addSysUnc,
       createSysDt = createSysDt,
       # necessary
       getErrTypes = getErrTypes,
       getSignature = getSignature,
       map = map)

}
