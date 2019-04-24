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


  createSysDt <- function(expIds, vals, uncs) {

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
    if (! dataref %in% names(expDt))
      stop(paste0("Column ", dataref, " specified as reference data ",
                  "must be present in expDt even if all ERRTYPE == sys-abs (in this case reference data is not used)"))
    setkey(assocDt, cat, lev)
    extSysDt <- merge(sysDt[, list(EXPID = EXPID, IDX2 = IDX, ERRTYPE = ERRTYPE)],
                      assocDt[, list(EXPID = expid, cat, lev)], by = "EXPID")
    uniqueCats <- unique(extSysDt[, cat])
    extExpDt <- expDt[, c("IDX", "EXPID", "DIDX", dataref, uniqueCats), with=FALSE]

    setnames(extExpDt, "IDX", "IDX1")
    extExpDt <- melt(extExpDt, id.vars = c("IDX1", "EXPID", "DIDX", dataref),
                     measure.vars = uniqueCats, variable.name = "cat", value.name = "lev")

    resDt <- merge(extExpDt[, list(IDX1, cat, lev, REFDATA = get(dataref))],
                   extSysDt[, list(IDX2, cat, lev, ERRTYPE)], by=c("cat", "lev"))
    resDt[, cat := NULL]
    resDt[, lev := NULL]
    resDt[, X := ifelse(ERRTYPE == "sys-abs", 1, REFDATA)]
    resDt[, ERRTYPE := NULL]
    resDt[, REFDATA := NULL]

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
