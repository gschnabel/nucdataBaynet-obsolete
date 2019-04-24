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


#' Verify Structure of SysDt
#'
#' @param sysDt TODO
#' @param before.prep TODO
#' @param print.info TODO
#'
#' @return TODO
#' @export
#'
isValidSysDt <- function(sysDt, range.check = TRUE, print.info = TRUE) {

  checkFun <- function(...) {
    checkExprs <- sapply(substitute(list(...)), deparse)[-1]
    passed <- sapply(list(...), isTRUE)
    if (!all(passed)) {
      cat(paste0("Failed checks: ",
                 paste0(checkExprs[!passed], collapse=", ")), "\n")
    }
    all(passed)
  }

  checkFun(is.data.table(sysDt),
           "EXPID" %in% colnames(sysDt),
           "DIDX" %in% colnames(sysDt),
           "GPTYPE" %in% colnames(sysDt),
           "ERRTYPE" %in% colnames(sysDt),
           "DATA" %in% colnames(sysDt),
           "UNC" %in% colnames(sysDt),
           is.character(sysDt$EXPID),
           is.integer(sysDt$DIDX),
           is.character(sysDt$ERRTYPE),
           is.numeric(sysDt$DATA),
           is.numeric(sysDt$UNC)) &&
    (!range.check ||
       checkFun("IDX" %in% colnames(sysDt),
                is.integer(sysDt$IDX),
                min(sysDt$IDX) == 1,
                max(sysDt$IDX) == nrow(sysDt)))
}


#' @name getDt_setDt_x
#' @rdname getDt_setDt_x
#'
#' @title Convenience Functions to Set and Get Columns in a Datatable
#'
#' @param dt a datatable
#' @param x a vector with new values of the column
#' @param idx if \code{NULL}, \code{x} will replace the complete column, otherwise
#' only the values associated with the indices in \code{idx} will be changed.
#'
#' @note Datatables are modified by reference, i.e. the datatable \code{dt} passed as argument is modified!
#'
NULL

#' @rdname getDt_setDt_x
#' @export
setDt_DATA <- function(dt, x, idx = NULL) {

  if (is.null(idx))
    idx <- seq_len(nrow(dt))
  setkey(dt, IDX)
  dt[J(idx), DATA := x]
}


#' @rdname getDt_setDt_x
#' @export
getDt_DATA <- function(dt, idx = NULL) {

  if (is.null(idx))
    idx <- seq_len(nrow(dt))
  setkey(dt, IDX)
  dt[J(idx), DATA]
}


#' @rdname getDt_setDt_x
#' @export
getDt_REFDATA <- function(dt, idx = NULL) {

  if (is.null(idx))
    idx <- seq_len(nrow(dt))
  setkey(dt, IDX)
  dt[J(idx), REFDATA]
}


#' @rdname getDt_setDt_x
#' @export
setDt_UNC <- function(dt, x, idx = NULL) {

  if (is.null(idx))
    idx <- seq_len(nrow(dt))
  setkey(dt, IDX)
  dt[J(idx), UNC := x]
}


#' @rdname getDt_setDt_x
#' @export
getDt_UNC <- function(dt, idx = NULL) {

  if (is.null(idx))
    idx <- seq_len(nrow(dt))
  setkey(dt, IDX)
  dt[J(idx), UNC]
}


