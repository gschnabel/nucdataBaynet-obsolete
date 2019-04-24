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


#' Construct Cholesky Decomposition
#'
#' @param x TODO
#'
#' @return TODO
#' @export
#'
makeCholesky <- function(x) {

  Cholesky(x, LDL = FALSE, perm = TRUE)
}


#' Title
#'
#' @param D TODO
#' @param S TODO
#' @param P TODO
#'
#' @return TODO
#' @export
#'
makeCholZ <- function(D, S, P) {

  Z <- forceSymmetric(crossprod(S, solve(D, S)) + solve(P))
  makeCholesky(Z)
}


#' Multiply an Inverse Covariance Matrix From Right
#'
#' @param D TODO
#' @param S TODO
#' @param P TODO
#' @param x TODO
#'
#' @return TODO
#' @import Matrix
#' @export
#'
mult_invCov_x <- function(x, D, S, P = NULL, cholZ = NULL) {

  stopifnot(grepl("Matrix$", class(D)),
            is.null(P) || grepl("Matrix$", class(P)),
            grepl("Matrix$", class(S)))

  invDx <- solve(D, x)
  if (is.null(cholZ)) cholZ <- makeCholZ(D, S, P)
  invDx - solve(D, S %*% solve(cholZ, crossprod(S, invDx)))
}


#' Multiply an Inverse Covariance Matrix from Both Sides
#'
#' @param D TODO
#' @param S TODO
#' @param P TODO
#' @param x TODO
#'
#' @return TODO
#' @import Matrix
#' @export
#'
mult_xt_invCov_x <- function(x, D, S, P = NULL, cholZ = NULL) {

  stopifnot(grepl("Matrix$", class(D)),
            is.null(P) || grepl("Matrix$", class(P)),
            grepl("Matrix$", class(S)))

  invDx <- solve(D, x)
  if (is.null(cholZ)) cholZ <- makeCholZ(D, S, P)
  crossprod(x, (invDx - solve(D, S %*% solve(cholZ, crossprod(S, invDx)))))
}


#' Trace of Inverse Covariance Matrix multplied with derivative of P
#'
#' @param dP TODO
#' @param D TODO
#' @param S TODO
#' @param P TODO
#' @param cholZ TODO
#'
#' @return TODO
#' @export
#'
tr_invCov_SdPtS <- function(dP, D, S, P = NULL, cholZ = NULL) {

  invDS <- solve(D, S)
  if (is.null(cholZ)) cholZ <- makeCholZ(D, S, P)
  pMat <- as(cholZ, "pMatrix")
  tSinvDS <- crossprod(S, invDS)
  invLtSinvDS <- solve(cholZ, pMat %*% tSinvDS, system = "L")
  print(rowSums(t(invLtSinvDS^2) * dP)) # debug
  colSums(S * invDS) * dP - rowSums(t(invLtSinvDS^2) * dP)
}



#' Calculate Chisquare
#'
#' @param D TODO
#' @param S TODO
#' @param P TODO
#' @param expdata TODO
#' @param moddata TODO
#'
#' @return TODO
#' @export
#'
chisquare <- function(d, D, S, P = NULL, cholZ = NULL) {

  as.vector(mult_xt_invCov_x(d, D, S, P, cholZ))
}


#' Logarithm of Determinant of Covariance matrix
#'
#' @param D TODO
#' @param S TODO
#' @param P TODO
#'
#' @return TODO
#' @export
#'
logDetCov <- function(D, S, P, cholZ = NULL) {

  stopifnot(grepl("Matrix$", class(D)),
            grepl("Matrix$", class(P)),
            grepl("Matrix$", class(S)))
  tmp <- determinant(D)
  stopifnot(tmp$sign == 1)
  logDetD <- tmp$modulus
  tmp <- determinant(P)
  stopifnot(tmp$sign == 1)
  logDetP <- tmp$modulus

  if (is.null(cholZ)) {
    Z <- forceSymmetric(crossprod(S, solve(D, S)) + solve(P))
    tmp <- determinant(crossprod(S, solve(D, S)) + solve(P))
  } else {
    tmp <- determinant(cholZ)
    tmp$modulus <- 2 * tmp$modulus
  }
  stopifnot(tmp$sign == 1)
  logDetZ <- tmp$modulus

  as.vector(logDetP + logDetD + logDetZ)
}


#' Calculate Log-Likelihood
#'
#' @param D TODO
#' @param S TODO
#' @param P TODO
#' @param expdata TODO
#' @param moddata TODO
#'
#' @return TODO
#' @export
#'
logLike <- function(d, D, S, P, cholZ = NULL) {

  if (is.null(cholZ)) cholZ <- makeCholZ(D, S, P)
  chisqr <- chisquare(d, D, S, P, cholZ)
  logDet <- logDetCov(D, S, P, cholZ)
  -0.5 * (chisqr + logDet + length(d) * log(2*pi))
}



#' Calculate Gradient of Log-Likelihood
#'
#' @param D TODO
#' @param S TODO
#' @param expdata TODO
#' @param moddata  TODO
#' @param sysDt  TODO
#' @param idxDt  TODO
#' @param sysCompHandler TODO
#'
#' @return TODO
#' @export
#'
gradLogLike <- function(d, D, S, cholZ = NULL, sysDt, gpDt = NULL, statUncIdx = NULL,
                                  sysUncIdx = NULL, gpIdx = NULL, sysCompHandler) {

  stopifnot("ddiMatrix" %in% class(D))

  if (is.null(statUncIdx)) statUncIdx <- seq_len(nrow(D))
  if (is.null(sysUncIdx)) sysUncIdx <- seq_len(nrow(sysDt))
  if (!is.null(gpDt) && is.null(gpIdx)) gpIdx <- seq_len(nrow(gpDt))

  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)

  # derivative of diagonal uncertainties
  gradChisqr_wrtStatUnc <- NULL
  gradLogDet_wrtStatUnc <- NULL
  gradLogLike_wrtStatUnc <- NULL

  gradChisqr_wrtSysUnc <- NULL
  gradLogDet_wrtSysUnc <- NULL
  gradLogLike_wrtSysUnc <- NULL

  gradChisqr_wrtHyp <- NULL
  gradLogDet_wrtHyp <- NULL
  gradLogLike_wrtHyp <- NULL

  setkey(sysDt, IDX)

  # helpful quantitiesn
  if (is.null(cholZ)) cholZ <- makeCholZ(D, S, P)
  # for chisquare derivative
  z1 <- mult_invCov_x(d, D, S, P, cholZ)
  z2 <- crossprod(S, z1)
  # for logDetCov derivative
  cholD <- makeCholesky(D)
  pMat <- as(cholD, "pMatrix")
  pMatS <- pMat %*% S
  if (is.null(dim(pMatS)))
    pMatS <- t(as(pMatS, "sparseMatrix"))

  L1 <- t(solve(cholD, pMatS, system = "L"))
  invDS <- solve(cholD, S)
  pMat <- as(cholZ, "pMatrix")
  tSinvDS <- crossprod(S, invDS)
  L2 <- t(solve(cholZ, pMat %*% tSinvDS, system = "L"))

  if (length(sysUncIdx) > 0) {
    gradChisqr_wrtSysUnc <- (-1) * sysCompHandler$grad_xtCovx_wrtUnc(z2, sysDt, idx = sysUncIdx, ret.mat=TRUE)
    gradLogDet_wrtSysUnc <- sysCompHandler$tr_LtdCovL_wrtUnc(L1, L2, sysDt, idx = sysUncIdx, ret.mat = TRUE)
    gradLogLike_wrtSysUnc <- (-0.5) * (gradChisqr_wrtSysUnc + gradLogDet_wrtSysUnc)
  }

  if (!is.null(gpDt) && length(gpIdx) > 0) {
    # gradChisqr_wrtHyp <- (-1) * sysCompHandler$grad_xtCovx_wrtHyp(z, sysDt, gpDt, idx = idx, ret.mat=TRUE)
    # gradLogDet_wrtHyp <- sysCompHandler$tr_LtdCovL_wrtHyp(L1, L2, sysDt, gpDt, idx = idx, ret.mat = TRUE)
    gradLogLike_wrtHyp <- sysCompHandler$gradLogLike_wrtHyp(z2, L1, L2, sysDt, gpDt, idx = gpIdx, ret.mat = TRUE)
  }

  # derivative wrt statistical uncertainties
  if (length(statUncIdx) > 0) {
    L3 <- solve(cholZ, pMat %*% t(invDS), system = "L")  # necessary for derivative wrt D aka statistical uncertainties
    invD <- solve(D)
    deriveD <- 2*sqrt(diag(D))
    gradChisqr_wrtStatUnc <- (-1) * deriveD * as.vector(z1)^2
    gradLogDet_wrtStatUnc <- diag(invD) * deriveD - colSums(L3 * t(deriveD * t(L3)))
    gradLogLike_wrtStatUnc <- (-0.5) * (gradChisqr_wrtStatUnc + gradLogDet_wrtStatUnc)
  }

  list(gradChisqr_wrtSysUnc = gradChisqr_wrtSysUnc,
       gradLogDet_wrtSysUnc = gradLogDet_wrtSysUnc,
       gradLogLike_wrtSysUnc = gradLogLike_wrtSysUnc,
       # gradChisqr_wrtHyp = gradChisqr_wrtHyp,
       # gradLogDet_wrtHyp = gradLogDet_wrtHyp,
       # gradLogLike_wrtHyp = (-0.5) * (gradChisqr_wrtHyp + gradLogDet_wrtHyp)
       gradLogLike_wrtHyp = gradLogLike_wrtHyp,
       gradLogLike_wrtStatUnc = gradLogLike_wrtStatUnc
  )
}



#' Create Functions for ML Optimization
#'
#' @return List with functions \code{logLike} and \code{gradLogLike}
#' @export
#'
createMLOptimFuns <- function() {

  # private variables
  this <- list(

    sysCompHandler = NULL,

    d = NULL,
    D = NULL,
    S = NULL,
    P = NULL,

    expDt = NULL,
    sysDt = NULL,
    gpDt = NULL,

    statUncIdx = NULL,
    sysUncIdx = NULL,
    gpIdx = NULL,

    idxSetStat = NULL,
    idxSetSys = NULL,
    idxSetGp = NULL
  )

  # setter/getter functions

  setDts <- function(expDt = NULL, sysDt = NULL, gpDt = NULL, sysCompHandler) {
    stopifnot(all(c("IDX","ADJUSTABLE","DATA","REFDATA") %in% names(expDt)),
              all(c("IDX","ADJUSTABLE","DATA","REFDATA") %in% names(sysDt)),
              is.null(gpDt) || all(c("IDX","ADJUSTABLE") %in% names(gpDt)))
    this$expDt <<- copy(expDt)
    this$sysDt <<- copy(sysDt)
    this$gpDt <<- copy(gpDt)
    setkey(this$expDt, IDX)
    setkey(this$sysDt, IDX)
    if (!is.null(this$gpDt))
      setkey(this$gpDt, IDX)
    this$statUncIdx <<- which(this$expDt$ADJUSTABLE)
    this$sysUncIdx <<- which(this$sysDt$ADJUSTABLE)
    if (!is.null(this$gpDt))
      this$gpIdx <<- which(this$gpDt$ADJUSTABLE)

    this$idxSetStat <<- seq_along(this$statUncIdx)
    this$idxSetSys <<- length(this$idxSetStat) + seq_along(this$sysUncIdx)
    this$idxSetGp <<- length(this$idxSetStat) + length(this$idxSetSys) + seq_along(this$gpIdx)

    this$S <<- sysCompHandler$map(this$expDt, this$sysDt, ret.mat = TRUE)
    this$P <<- sysCompHandler$cov(this$sysDt, this$gpDt, ret.mat = TRUE)
    this$d <<- getDt_DATA(this$expDt) - getDt_REFDATA(this$expDt) -
      this$S %*% (getDt_DATA(this$sysDt) - getDt_REFDATA(this$sysDt))

    this$sysCompHandler <<- sysCompHandler
  }

  # functions

  updateDts <- function(x) {

    stopifnot(length(x) == length(this$idxSetStat) + length(this$idxSetSys) + length(this$idxSetGp))
    setkey(this$expDt, IDX)
    this$expDt[J(this$statUncIdx), UNC := x[this$idxSetStat]]
    setkey(this$sysDt, IDX)
    this$sysDt[J(this$sysUncIdx), UNC := x[this$idxSetSys]]
    if (!is.null(this$gpDt)) {
      setkey(this$gpDt, IDX)
      this$gpDt[J(this$gpIdx), PARVAL := x[this$idxSetGp]]
    }
    this$P <<- this$sysCompHandler$cov(this$sysDt, this$gpDt, ret.mat = TRUE)
    this$D <<- Diagonal(x = getDt_UNC(this$expDt)^2)
  }

  getModifiedDts <- function(x) {

    updateDts(x)
    list(statUncIdx = this$statUncIdx,
         sysUncIdx = this$sysUncIdx,
         gpIdx = this$gpIdx,
         expDt = copy(this$expDt),
         sysDt = copy(this$sysDt),
         gpDt = copy(this$gpDt))
  }

  locLogLike <- function(x) {

    updateDts(x)
    logLike(this$d, this$D, this$S, this$P)
  }

  locGradLogLike <- function(x) {

    updateDts(x)
    gradStruc <- gradLogLike(this$d, this$D, this$S,
                             sysDt = this$sysDt,
                             gpDt = this$gpDt,
                             statUncIdx = this$statUncIdx,
                             sysUncIdx = this$sysUncIdx,
                             gpIdx = this$gpIdx,
                             sysCompHandler = this$sysCompHandler)
    with(gradStruc, c(
      gradLogLike_wrtStatUnc[this$statUncIdx],
      gradLogLike_wrtSysUnc[this$sysUncIdx],
      gradLogLike_wrtHyp[this$gpIdx]
    ))
  }

  list(setDts = setDts,
       getModifiedDts = getModifiedDts,
       # computation functions
       logLike = locLogLike,
       gradLogLike = locGradLogLike)
}







