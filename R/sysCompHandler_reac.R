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


#' Create Systematic Component Handler for Reactions
#'
#' @param subentList TODO
#'
#' @return TODO
#' @import data.table Matrix
#' @export
#'
createSysCompReacHandler <- function(subentList) {

  # initialization
  assocDt <- rbindlist(lapply(subentList, function(x) {
    curReac <- x$BIB$REACTION
    if (is.list(curReac) || length(curReac)!=1) NULL
    else if (is.null(x$DATA$TABLE$EN)) NULL else {
      data.table(EXPID = x$ID, DIDX = seq_len(nrow(x$DATA$TABLE)),
                 EN = x$DATA$TABLE$EN, REAC = reacStrucToStr(parseReacExpr(curReac)))
    }}))
  setkey(assocDt, REAC, EXPID, DIDX)

  uniqueReacs <- unique(assocDt$REAC)
  mapAssignmentDt <- data.table(IDX = integer(0),
                                EXPID = character(0),
                                REAC = character(0),
                                VALS = list(),
                                UNCS = list(),
                                HNDNAME = character(0),
                                OPTS = list())

  mapList <- list()  # list of lists with mapper functions

  # administration

  addMap <- function(hndname, mapper) {
    stopifnot(is.function(mapper))
    mapList[[hndname]] <<- mapper
  }

  getMapList <- function() {
    mapList
  }

  assignMapToReac <- function(hndname, reac, vals, uncs, opts = list()) {

    stopifnot(hndname %in% names(mapList))
    stopifnot(length(vals) == length(uncs))
    stopifnot(is.list(opts))

    reacStruc <- parseReacExpr(reac)
    stopifnot(!is.null(reacStruc))
    stopifnot(isTRUE(reacStruc$quantspec == ",SIG"))

    reacDescStr <- with(reacStruc, paste0(projectile, ",", process))
    newSysExpIdStem <- sprintf("REACEXP-%s-", reacDescStr)
    reacIdx <- mapAssignmentDt[, sum(grepl(newSysExpIdStem, EXPID, fixed = TRUE)) + 1]

    newIdx <- max(c(mapAssignmentDt$IDX,0)) + 1
    newRowDt <- data.table(IDX = newIdx,
                           EXPID = sprintf("REACEXP-%s-%02d",
                                           reacDescStr, reacIdx),
                           REAC = reac,
                           VALS = list(vals),
                           UNCS = list(uncs),
                           HNDNAME = hndname, OPTS = list(opts))
    mapAssignmentDt <<- rbind(mapAssignmentDt, newRowDt)
  }

  getMapAssignment <- function() {
    mapAssignmentDt
  }

  createSysDt <- function() {
    resDt <- mapAssignmentDt[, {
      list(DIDX = seq_along(VALS[[1]]),
           ERRTYPE = HNDNAME,
           GPTYPE = NA_character_,
           DATA = VALS[[1]],
           UNC = UNCS[[1]]
      )
    }, by="EXPID"]
    resDt
  }


  getSignature <- function() {
    "reactionMapHandler"
  }


  getErrTypes <- function() {
    names(mapList)
  }

  # mapping

  map <- function(expDt, sysDt, ret.mat) {

    # assumption: expDt contains only one virtual experiment
    setkey(mapAssignmentDt, EXPID)

    extExpDt <- merge(assocDt, expDt[,list(IDX, EXPID, DIDX)],
                      by = c("EXPID", "DIDX"))
    setkey(extExpDt, REAC)

    resDt <- sysDt[, {
      curMapAssignDt <- mapAssignmentDt[J(.BY[["EXPID"]])]
      stopifnot(nrow(curMapAssignDt) == 1)
      curExpDt <- extExpDt[J(curMapAssignDt$REAC), list(IDX, EN)]
      curMapper <- mapList[[curMapAssignDt$HNDNAME]]
      stopifnot(is.function(curMapper))
      if (!is.na(curExpDt$EN[1])) {
        curResDt <- curMapper(curExpDt$EN,
                              opts = curMapAssignDt$OPTS[[1]],
                              ret.mat = FALSE)
        globalIdx1 <- curExpDt$IDX
        globalIdx2 <- IDX
        curResDt[, IDX1 := globalIdx1[IDX1]]
        curResDt[, IDX2 := globalIdx2[IDX2]]
        curResDt
      } else {
        data.table(IDX1 = integer(0), IDX2 = integer(0), X = numeric(0))
      }
    }, by="EXPID"]
    resDt[, EXPID := NULL]
    resDt <- resDt[X != 0, ]
    resDims <- c(nrow(expDt), nrow(sysDt))
    returnSparseMatrix(resDt, resDims, ret.mat)
  }


  list(addMap = addMap,
       getMapList = getMapList,
       assignMapToReac = assignMapToReac,
       getMapAssignment = getMapAssignment,
       createSysDt = createSysDt,
       # must have properties
       getSignature = getSignature,
       getErrTypes = getErrTypes,
       map = map)
}
