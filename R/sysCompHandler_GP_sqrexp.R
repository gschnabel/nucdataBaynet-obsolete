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


#' Create Handler for Systematic GP Components
#'
#' @return TODO
#' @import data.table Matrix
#' @export
#'
createSysCompGPHandler <- function() {

  gpDt <- NULL

  addGP <- function(expid, sigma, len, nugget) {

    newDt <- data.table(EXPID = expid,
                        GPTYPE = "sqrexp",
                        PARNAME = c("sigma", "len", "nugget"),
                        PARVAL = c(sigma, len, nugget))
    proposedDt <- rbind(gpDt, newDt)
    stopifnot(anyDuplicated(proposedDt) == 0)
    gpDt <<- proposedDt
  }

  createGPDt <- function() {
    gpDt
  }

  updateSysDt <- function(sysDt) {

    uniqDt <- unique(gpDt[, list(EXPID, GPTYPE)])
    setkey(uniqDt, "EXPID")
    sysDt[, GPTYPE := uniqDt[J(sysDt$EXPID), GPTYPE]]
    sysDt[]
  }

  getSignature <- function() {
    "GPhandler"
  }

  getGPTypes <- function() {
    "sqrexp"
  }

  cov <- function(sysDt, gpDt, ret.mat) {

    setkey(gpDt, EXPID)
    stopifnot(is.data.table(sysDt),
              is.data.table(gpDt))
    stopifnot("IDX" %in% names(sysDt))

    resDt <- sysDt[GPTYPE == "sqrexp", {

      curGPSpec <- gpDt[J(.BY[["EXPID"]])]
      curSigma <- curGPSpec[PARNAME == "sigma", PARVAL]
      curLen <- curGPSpec[PARNAME == "len", PARVAL]
      curNugget <- curGPSpec[PARNAME == "nugget", PARVAL]

      covMat <- curSigma^2 * exp(-0.5/curLen^2 * (outer(EN, EN, `-`))^2)
      diag(covMat) <- diag(covMat) + curNugget^2
      list(IDX1 = IDX, IDX2 = rep(IDX, each=length(IDX)), X = as.vector(covMat))
    }, by="EXPID"]
    resDt[,EXPID := NULL]
    resDims <- c(nrow(sysDt), nrow(sysDt))
    returnSparseMatrix(resDt, resDims, ret.mat)
  }


  basicGradfun_wrtHyp <- function(x = NULL, L1 = NULL, L2 = NULL,
                                  sysDt, gpDt = NULL, idx = NULL, ret.mat, what) {

    stopifnot("IDX" %in% names(gpDt))
    if (is.null(idx)) idx <- seq_len(nrow(gpDt))

    setkey(gpDt, IDX)
    selGpDt <- gpDt[idx,]

    setkey(gpDt, EXPID, PARNAME)
    setkey(sysDt, EXPID)

    # precomputing L1 and L2
    # if (what == "tr_LtdCovL_wrtHyp") {
    #   L1tL1 <- forceSymmetric(L1 %*% t(L1))
    #   L2tL2 <- forceSymmetric(L2 %*% t(L2))
    # }

    gradVec <- rep(0, nrow(gpDt))
    for (i in seq_len(nrow(selGpDt))) {

      curGpIdx <- selGpDt[i, IDX]
      curExpId <- selGpDt[i, EXPID]
      idxEnDt <- sysDt[J(curExpId), list(IDX, EN)]

      tmp <- gpDt[J(curExpId, c("sigma", "len", "nugget")), PARVAL]
      curSigma <- tmp[1]
      curLen <- tmp[2]
      curNugget <- tmp[3]

      curParname <- selGpDt[i, PARNAME]
      curParval <- selGpDt[i, PARVAL]

      if (curParname == "sigma") {
        z <- (-0.5) * (outer(idxEnDt$EN, idxEnDt$EN, `-`))^2
        dCovMat <- 2*curSigma * exp(z/curLen^2)
      } else if (curParname == "len") {
        z <- (-0.5) * (outer(idxEnDt$EN, idxEnDt$EN, `-`))^2
        dCovMat <- curSigma^2 * exp(z/curLen^2) * ((-2)*z)/curLen^3
      } else {
        dCovMat <- matrix(0, nrow = length(idxEnDt$EN), ncol = length(idxEnDt$EN))
      }

      gradVec[curGpIdx] <- switch(what,
                                  grad_xtCovx_wrtHyp = {
                                    redx <- x[idxEnDt$IDX]
                                    as.vector(t(redx) %*% dCovMat %*% redx)
                                  },
                                  tr_LtdCovL_wrtHyp = {
                                    redL1 <- as.matrix(L1[idxEnDt$IDX,,drop=FALSE])  # as.matrix makes it slightly faster
                                    redL2 <- as.matrix(L2[idxEnDt$IDX,,drop=FALSE])  # at least for matrices of size 16x16
                                    sum(redL1 * (dCovMat %*% redL1)) - sum(redL2 * (dCovMat %*% redL2))
                                    # another approach
                                    # redL1tL1 <- L1tL1[idxEnDt$IDX,idxEnDt$IDX,drop=FALSE] # seems to be slower
                                    # redL2tL2 <- L2tL2[idxEnDt$IDX,idxEnDt$IDX,drop=FALSE]
                                    # sum(redL1tL1 * dCovMat) - sum(redL2tL2 * dCovMat)
                                  },
                                  gradLogLike_wrtHyp = {
                                    redx <- x[idxEnDt$IDX]
                                    redL1 <- as.matrix(L1[idxEnDt$IDX,,drop=FALSE])
                                    redL2 <- as.matrix(L2[idxEnDt$IDX,,drop=FALSE])
                                    (-0.5) * (sum(redL1 * (dCovMat %*% redL1)) - sum(redL2 * (dCovMat %*% redL2)) -
                                                as.vector(t(redx) %*% dCovMat %*% redx))
                                  })
    }
    resDt <- data.table(IDX1 = idx, IDX2 = 1L, X = gradVec[idx])

    resDims <- c(1, nrow(gpDt))
    returnSparseMatrix(resDt, resDims, ret.mat)
  }


  grad_xtCovx_wrtHyp <- function(x, sysDt, gpDt, idx = NULL, ret.mat) {

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

  list(addGP = addGP,
       createGPDt = createGPDt,
       updateSysDt = updateSysDt,
       getSignature = getSignature,
       getGPTypes = getGPTypes,
       cov = cov,
       grad_xtCovx_wrtHyp = grad_xtCovx_wrtHyp,
       tr_LtdCovL_wrtHyp = tr_LtdCovL_wrtHyp,
       gradLogLike_wrtHyp = gradLogLike_wrtHyp)
}
