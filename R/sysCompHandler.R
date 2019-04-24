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


#' Create Systematic Component Handler
#'
#' Creates a handler (i.e. list of functions) to operate with systematic components.
#'
#' @return A list with functions to manage systematic components, see \code{Details}
#' @details
#' The two main responsibilities of the handler are (1) to construct a mapping matrix \code{S}
#' that maps the systematic components to the measurements, and (2) to construct a covariance matrix
#' \code{P} associated with the systematic components.
#' The respective functions are
#' \describe{
#'   \item{\code{map(expDt, sysDt, ret.mat)}}{
#'     The datatable \code{expDt} must contain the columns \code{IDX, EXPID, DIDX, DATA, UNC}.
#'     The datatable \code{sysDt} must contain the columns \code{IDX, EXPID, DIDX, ERRTYPE, GPTYPE, DATA, UNC}.
#'     If the logical flag \code{ret.mat} is \code{TRUE}, a sparse matrix will be returned,
#'     otherwise a datatable with columns \code{IDX1, IDX2, X}.
#'   }
#'   \item{\code{cov(sysDt, gpDt = NULL, ret.mat)}}{
#'     The datatable \code{sysDt} is as above.
#'     The datatable \code{gpDt} must contain the columns \code{IDX, EXPID, GPTYPE, PARNAME, PARVAL}.
#'     The logical flag \code{ret.mat} is as above.
#'     The datatable \code{gpDt} is only needed if \code{sysDt} contains non-\code{NA} elements in \code{GPTYPE}.
#'   }
#' }
#' The function \code{map} is dependent on handlers that know how to map specific error types
#' given in the column \code{ERRTYPE} of datatable \code{sysDt}.
#' The following functions are available to manage handlers
#' \describe{
#'   \item{\code{addHandler(handler)}}{
#'     The list \code{handler} must contain the functions \code{getSignature, getErrTypes, map}.
#'   }
#'   \item{\code{getHandlerList()}}{
#'     Returns a list of the handlers registered by \code{addHandler}.
#'   }
#'   \item{\code{getHandlerNames()}}{
#'     Returns a character vector with the names of the handlers.
#'   }
#' }
#' The function \code{cov} by itself is able to provide the elements of the covariance matrix associated with rows in \code{sysDt} that
#' have a \code{NA} value in column \code{GPTYPE}.
#' For the blocks of the covariance matrix associated with rows with a non-\code{NA} as \code{GPTYPE}, appropriate subhandlers must be
#' registered that know how to calculate the covariance blocks for the specific Gaussian process types.
#' The following functions are available to manage GP handlers:
#' \describe{
#'   \item{\code{addGPHandler(handler)}}{
#'     The list \code{handler} must contain the functions \code{getSignature, getGPTypes, cov}.
#'   }
#'   \item{\code{getGPHandlerList()}}{
#'     Returns a list of the GP handlers registered by \code{addGPHandler}
#'   }
#'   \item{\code{getGPHandlerNames()}}{
#'     Returns a character vector with the names of the GP handlers.
#'   }
#' }
#' Maximum likelihood optimization for large scale problems requires some mathematical tricks to be feasible.
#' The following functions are helpful for the task:
#' \describe{
#'   \item{\code{grad_xtCovx_wrtUnc(x, sysDt, idx = NULL, ret.mat)}}{
#'     Evaluates the gradient of the expression \code{t(x) \%*\% K \%*\% x} with respect to partial derivatives of the
#'     uncertainties whose squares constitute the diagonal elements of the covariance matrix K.
#'     The matrix K in this expression is as returned by function \code{cov} explained above.
#'     The integer vector \code{idx} allows to calculate the gradient with respect to a subset of the uncertainties,
#'     other elements in the resulting gradient vector are set to zero.
#'     The meaning of \code{ret.mat} is as above.
#'   }
#'   \item{\code{grad_xtCovx_wrtHyp(x, sysDt, gpDt, idx = NULL, ret.mat)}}{
#'     Evaluates the gradient of the expression \code{t(x) \%*\% K \%*\% x} with respect to partial derivatives of the
#'     hyperparameters defined in \code{gpDt}.
#'     The integer vector \code{idx} allows to calculate the gradient with respect to a subset of the hyperparameters,
#'     other elements in the resulting gradient vector are set to zero.
#'     The meaning of \code{ret.mat} is as above.
#'   }
#'   \item{\code{tr_LtCovL_wrtUnc(L1, L2, sysDt, ret.mat)}}{
#'     Given matrices L1 and L2, evaluates the expression
#'     \eqn{\textrm{tr}(L_1^T \partial K/\partial{\delta_i} L_1) + L_2^T \partial K / \partial \delta_i L_2}
#'     where \eqn{\partial K / \partial \delta_i} is the derivative with respect to an uncertainty \eqn{\delta_i}.
#'     The function returns a vector where the \eqn{i^\textrm{th}} element corresponds to the derivative wrt \eqn{\delta_i}.
#'   }
#' }
#'
#'
#' @import data.table Matrix
#' @export
#'
#'
createSysCompHandler <- function() {

  handlerList <- list()
  handlerErrTypeAssignDt <- NULL

  gpHandlerList <- list()
  handlerGPTypeAssignDt <- NULL

  # named list, name given by signature
  # each list element contains functions map, getSignature, getErrTypes
  # should also include: covDeriv, SdPSd, etc.

  addHandler <- function(handler) {
    stopifnot(is.list(handler))
    stopifnot(c("map", "getErrTypes",
                "getSignature") %in% names(handler))
    with(handler, stopifnot(is.function(map), is.function(getSignature),
                            is.function(getErrTypes)))
    newHandlerList <- c(handlerList, list(handler))
    newNames <- c(names(handlerList), handler$getSignature())
    names(newHandlerList) <- newNames
    handlerList <<- newHandlerList

    handlerErrTypeAssignDt <<- rbindlist(
      lapply(handlerList, function(x) {
        data.table(HNDNAME = x$getSignature(),
                   ERRTYPE = x$getErrTypes())
      }))
  }

  getHandlerList <- function() {
    handlerList
  }

  getHandlerNames <- function() {
    names(handlerList)
  }


  addGPHandler <- function(handler) {
    stopifnot(is.list(handler))
    stopifnot(c("cov", "getGPTypes",
                "getSignature") %in% names(handler))
    with(handler, stopifnot(is.function(map), is.function(getSignature),
                            is.function(getGPTypes)))
    newGPHandlerList <- c(gpHandlerList, list(handler))
    newNames <- c(names(gpHandlerList), handler$getSignature())
    names(newGPHandlerList) <- newNames
    gpHandlerList <<- newGPHandlerList

    handlerGPTypeAssignDt <<- rbindlist(
      lapply(gpHandlerList, function(x) {
        data.table(HNDNAME = x$getSignature(),
                   GPTYPE = x$getGPTypes())
      }))
  }

  getGPHandlerList <- function() {
    gpHandlerList
  }

  getGPHandlerNames <- function() {
    names(gpHandlerList)
  }

  createSysDt <- function() {

    sysDt <- rbindlist(lapply(handlerList, function(x) x$createSysDt()), fill = TRUE, idcol = "IDX")
    sysDt[, IDX := seq_len(.N)]
    sysDt[]
  }

  map <- function(expDt, sysDt, ret.mat) {

    # TODO: add isValidExpDt function
    stopifnot(isValidSysDt(sysDt, range.check = TRUE))
    extSysDt <- merge(sysDt, handlerErrTypeAssignDt, by = "ERRTYPE")
    if (nrow(extSysDt) == 0)
      stop("Either No systematic components in sysDt, or no handler available that knows how to deal with the existing ones")

    resDt <- extSysDt[,{
      curHandlerName <- HNDNAME[1]
      if (!all(HNDNAME == curHandlerName))
        stop("It is not allowed to have several handlers for the same virtual experiment (EXPID)")
      curHandler <- handlerList[[curHandlerName]]
      curRes <- curHandler$map(expDt, data.table(ERRTYPE, .SD), ret.mat = FALSE)
      if (nrow(curRes) > 0) curRes else list(IDX1 = integer(0), IDX2 = integer(0), X = numeric(0))
    }, by="ERRTYPE"]
    resDt[, ERRTYPE := NULL]
    resDt <- resDt[X != 0, ]

    resDims <- c(nrow(expDt), nrow(sysDt))
    returnSparseMatrix(resDt, resDims, ret.mat)
  }

  cov <- function(sysDt, gpDt = NULL, ret.mat = FALSE) {

    resDt <- sysDt[is.na(GPTYPE), list(IDX1 = IDX, IDX2 = IDX, X = UNC^2)]
    resDt <- resDt[X != 0, ]

    # if something is dependent on GP
    if (!is.null(gpDt)) {
      extSysDt <- merge(sysDt, handlerGPTypeAssignDt, by = "GPTYPE")
      gpResDt <- extSysDt[!is.na(GPTYPE), {
        curGPHandlerName <- HNDNAME[1]
        if (!all(HNDNAME == curGPHandlerName))
          stop("It is not allowed to have several handlers for the same virtual experiment (EXPID)")
        curGPHandler <- gpHandlerList[[curGPHandlerName]]
        curGPHandler$cov(data.table(GPTYPE, .SD), gpDt, ret.mat = FALSE)
      }, by = "GPTYPE"]
      gpResDt[, GPTYPE := NULL]
      resDt <- rbind(resDt, gpResDt)
    }

    resDims <- c(nrow(sysDt), nrow(sysDt))
    returnSparseMatrix(resDt, resDims, ret.mat)
  }

  # functions for derivatives

  grad_xtCovx_wrtUnc <- function(x, sysDt, idx = NULL, ret.mat) {

    setkey(sysDt, IDX)
    stopifnot(length(x) == nrow(sysDt))
    if (is.null(idx)) idx <- seq_along(x)
    # if a GPTYPE is present, the UNC column is not used
    idx <- idx[sysDt[idx, is.na(GPTYPE)]]
    resDt <- data.table(IDX1 = idx, IDX2 = 1L, X = 2*(x[idx]*x[idx]*sysDt$UNC[idx]))
    resDims <- c(nrow(sysDt),1)
    returnSparseMatrix(resDt, resDims, ret.mat)
  }


  tr_LtdCovL_wrtUnc <- function(L1, L2, sysDt, idx = NULL, ret.mat) {

    setkey(sysDt, IDX)
    if (is.null(idx)) idx <- seq_len(nrow(sysDt))
    # if a GPTYPE is present, the UNC column is not used
    idx <- idx[sysDt[idx, is.na(GPTYPE)]]
    dP <- 2*sysDt$UNC
    dP[-idx] <- 0

    resDt <- data.table(IDX1 = sysDt$IDX,
                        IDX2 = 1L,
                        X = rowSums(L1^2 * dP) - rowSums(L2^2 * dP))
    resDims <- c(nrow(sysDt), 1)
    returnSparseMatrix(resDt, resDims, ret.mat)
  }


  basicGradfun_wrtHyp <- function(x = NULL, L1 = NULL, L2 = NULL,
                                  sysDt, gpDt = NULL, idx = NULL, ret.mat, what) {

    if (is.null(idx)) idx <- seq_len(nrow(gpDt))
    extSysDt <- merge(sysDt, handlerGPTypeAssignDt, by = "GPTYPE")
    setkey(gpDt, IDX)
    extGpDt <- copy(gpDt)
    extGpDt[, deriveFlag := FALSE]
    extGpDt[J(idx), deriveFlag := TRUE]

    gpResDt <- extSysDt[!is.na(GPTYPE), {
      curGpType <- .BY[["GPTYPE"]]
      curGPHandlerName <- HNDNAME[1]
      if (!all(HNDNAME == curGPHandlerName))
        stop("It is not allowed to have several handlers for the same virtual experiment (EXPID)")
      curGPHandler <- gpHandlerList[[curGPHandlerName]]
      curGpDt <- extGpDt[GPTYPE == curGpType]
      curSysDt <- data.table(GPTYPE, .SD)

      switch(what,
             grad_xtCovx_wrtHyp = curGPHandler$grad_xtCovx_wrtHyp(x, curSysDt, curGpDt,
                                                                 curGpDt[deriveFlag==TRUE, IDX],
                                                                 ret.mat = FALSE),
             tr_LtdCovL_wrtHyp = curGPHandler$tr_LtdCovL_wrtHyp(L1, L2, curSysDt, curGpDt,
                                                                curGpDt[deriveFlag==TRUE, IDX],
                                                                ret.mat = FALSE),
             gradLogLike_wrtHyp = curGPHandler$gradLogLike_wrtHyp(x, L1, L2, curSysDt, curGpDt,
                                                                 curGpDt[deriveFlag==TRUE, IDX],
                                                                 ret.mat = FALSE) )
    }, by = "GPTYPE"]
    gpResDt[, GPTYPE := NULL]

    resDims <- c(nrow(gpDt),1)
    returnSparseMatrix(gpResDt, resDims, ret.mat)
  }


  grad_xtCovx_wrtHyp <- function(x, sysDt, gpDt, idx = NULL, ret.mat) {

    stopifnot(length(x) == nrow(sysDt))
    basicGradfun_wrtHyp(x = x, sysDt = sysDt, gpDt = gpDt,
                        idx = idx, ret.mat = ret.mat, what = "grad_xtCovx_wrtHyp")
  }


  tr_LtdCovL_wrtHyp <- function(L1, L2, sysDt, gpDt, idx = NULL, ret.mat) {

    basicGradfun_wrtHyp(L1 = L1, L2 = L2, sysDt = sysDt, gpDt = gpDt,
                        idx = idx, ret.mat = ret.mat, what = "tr_LtdCovL_wrtHyp")
  }


  gradLogLike_wrtHyp <- function(x, L1, L2, sysDt, gpDt, idx = NULL, ret.mat) {

    basicGradfun_wrtHyp(x = x, L1 = L1, L2 = L2, sysDt = sysDt, gpDt = gpDt,
                        idx = idx, ret.mat = ret.mat, what = "gradLogLike_wrtHyp")
  }


  list(# handler management
    addHandler = addHandler,
    getHandlerList = getHandlerList,
    getHandlerNames = getHandlerNames,
    addGPHandler = addGPHandler,
    getGPHandlerList = getGPHandlerList,
    getGPHandlerNames = getGPHandlerNames,
    createSysDt = createSysDt,
    # essential
    map = map, cov = cov,
    # helpers for derivative, mainly used by functions in marlikefuns.R
    grad_xtCovx_wrtUnc = grad_xtCovx_wrtUnc,
    tr_LtdCovL_wrtUnc = tr_LtdCovL_wrtUnc,
    grad_xtCovx_wrtHyp = grad_xtCovx_wrtHyp,
    tr_LtdCovL_wrtHyp = tr_LtdCovL_wrtHyp,
    gradLogLike_wrtHyp = gradLogLike_wrtHyp)
}

