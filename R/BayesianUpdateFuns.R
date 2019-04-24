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


#' Get Posterior of Systematic Components
#'
#' @param expDt datatable with experimental data
#' @param sysDt datatable with systematic components
#' @param gpDt datatable with Gaussian process components
#' @param sysCompHandler systematic component handler, see \link{createSysCompHandler}
#'
#' @return list with two elements \code{sysDt} and \code{cov} which contain
#' the posterior of the systematic components and the covariance matrix, respectively.
#' and covariance matrix
#' @export
#'
getPosteriorSys <- function(expDt, sysDt, gpDt = NULL, sysCompHandler) {

  setkey(expDt, IDX)
  S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)
  D <- Diagonal(x = getDt_UNC(expDt)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)

  invD <- solve(D)
  invP <- solve(forceSymmetric(P))
  p0 <-getDt_DATA(sysDt)
  expData <- getDt_DATA(expDt)
  refExpData <- getDt_REFDATA(expDt)
  refSysData <- getDt_REFDATA(sysDt)
  effData <- expData - refExpData + S %*% refSysData

  P1 <- solve(forceSymmetric(t(S) %*% invD %*% S) + invP)
  p1 <- as.vector(P1 %*% (invP %*% p0 + t(S) %*% invD %*% effData))

  newSysDt <- copy(sysDt)
  setDt_DATA(newSysDt, p1)
  setDt_UNC(newSysDt, sqrt(diag(P1)))
  list(sysDt = newSysDt, cov = P1)
}




#' Map Systematic Components to Experimental Grid
#'
#' @param expDt  datatable with experimental data
#' @param sysList list with components \code{sysDt} and potentially \code{cov}
#' @param sysCompHandler handler for systematic components as created by \link{createSysCompHandler}
#' @param ret.cov should the posterior matrix be computed?
#' @param idx vector of indices that specifies which systematic components should be mapped.
#' Default \code{NULL} means to map all of them.
#'
#' @return list with two components \code{expDt} and \code{cov}
#' @export
#'
mapSysToExp <- function(expDt, sysList, sysCompHandler, ret.cov = FALSE, idx = NULL) {

  if (is.null(idx)) idx <- seq_len(nrow(sysList$sysDt))
  sysDt <- sysList$sysDt
  S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)[,idx]
  refExpData <- getDt_REFDATA(expDt)
  refSysData <- getDt_REFDATA(sysDt, idx = idx)
  sysData <- getDt_DATA(sysDt, idx = idx)

  newExpDt <- copy(expDt)
  xsCovmat <- NULL
  setDt_DATA(newExpDt, as.vector(refExpData + S %*% (sysData - refSysData)))
  if (!is.null(sysList$cov)) {
    curCov <- sysList$cov[idx,idx]
    setDt_UNC(newExpDt, sqrt(rowSums((S %*% curCov) * S)))
    if (isTRUE(ret.cov)) {
      xsCovmat <- S %*% curCov %*% t(S)
    }
  }
  list(expDt = newExpDt[], cov = xsCovmat)
}




