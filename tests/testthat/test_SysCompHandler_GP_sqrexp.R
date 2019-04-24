context("testing GP functionality of sysCompHandler")

library(numDeriv)
library(testthat)
library(talysExforMapping)

subentHandler <- createSubentHandler(createDefaultSubentHandlerList())
exforHandler <- createExforHandler(subentHandler)

talysHandler <- createSysCompModelHandler()
talysHandler$setRef(exampleNeedsDt, exampleSparDt, exampleParamDt, exforHandler, exampleSubents)
talysHandler$setPrior(exampleParamDt)

sysDt <- talysHandler$createSysDt()
sysDt[, IDX := seq_len(.N)]

gpHandler <- createSysCompGPHandler()
sysDt[ERRTYPE == "talyspar_endep",
      gpHandler$addGP(.BY[["EXPID"]], 0.1, 5, 1e-3)
      , by="EXPID"]
gpDt <- gpHandler$createGPDt()
gpDt[, IDX := seq_len(.N)]
gpHandler$updateSysDt(sysDt)

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(talysHandler)
sysCompHandler$addGPHandler(gpHandler)

numericGradfunHyp <- function(x, idx, z) {

  setkey(sysDt, IDX)
  setkey(gpDt, IDX)
  cpGpDt <- copy(gpDt)
  cpGpDt[idx, PARVAL := PARVAL + x]
  curCov <- sysCompHandler$cov(sysDt, cpGpDt, ret.mat = TRUE)
  as.vector(t(z) %*% curCov %*% z)
}


numericCovfunHyp <- function(x, idx) {

  setkey(sysDt, IDX)
  setkey(gpDt, IDX)
  cpGpDt <- copy(gpDt)
  cpGpDt[idx, PARVAL := PARVAL + x]
  curCov <- sysCompHandler$cov(sysDt, cpGpDt, ret.mat = TRUE)
  as.matrix(curCov)
}


numericGradfunUnc <- function(x, idx, z) {

  setkey(sysDt, IDX)
  setkey(gpDt, IDX)
  cpSysDt <- copy(sysDt)
  cpSysDt[idx, UNC := UNC + x]
  curCov <- sysCompHandler$cov(cpSysDt, gpDt, ret.mat = TRUE)
  as.vector(t(z) %*% curCov %*% z)
}


numericGradfunLoglike <- function(x, idx, d, D, S, cholZ = NULL, parType) {

  setkey(sysDt, IDX)
  setkey(gpDt, IDX)
  cpSysDt <- copy(sysDt)
  cpGpDt <- copy(gpDt)

  if (parType == "sysUnc")
    cpSysDt[idx, UNC := UNC + x]
  else if (parType == "statUnc")
    D[idx,idx] <- (sqrt(D[idx,idx]) + x)^2
  else if (parType == "gpPar")
    cpGpDt[idx, PARVAL := PARVAL + x]
  else stop("unsupported parameter type")

  P <- sysCompHandler$cov(cpSysDt, cpGpDt, ret.mat = TRUE)
  logLike(d, D, S, P, cholZ)
}


test_that("grad_xtCovx_wrtHyp produces correct result", {

  d <- runif(nrow(sysDt))
  analyticGrad <- as.vector(sysCompHandler$grad_xtCovx_wrtHyp(d, sysDt, gpDt, ret.mat = TRUE))

  idcs <- c(1,2)
  numericGrad <- rep(0, nrow(gpDt))
  for (idx in idcs)
    numericGrad[idx] <- as.vector(jacobian(numericGradfunHyp, x = 0, idx = idx, z = d))
  expect_equal(analyticGrad[idcs], numericGrad[idcs])
  # subsetting
  analyticGrad <- as.vector(sysCompHandler$grad_xtCovx_wrtHyp(d, sysDt, gpDt, idx = idcs, ret.mat = TRUE))
  expect_equal(analyticGrad, numericGrad)
})


test_that("tr_LtdCovL produces correct result", {

  L1 <- matrix(runif(nrow(sysDt)*(nrow(sysDt)+5)), nrow(sysDt), nrow(sysDt) + 5)
  L2 <- matrix(runif(nrow(sysDt)*(nrow(sysDt)+5)), nrow(sysDt), nrow(sysDt) + 5)

  analyticGrad <- as.vector(sysCompHandler$tr_LtdCovL_wrtHyp(L1, L2, sysDt, gpDt, ret.mat = TRUE))
  idcs <- c(1,2)
  numericGrad <- rep(0, nrow(gpDt))
  for (idx in idcs) {
    dCov <- jacobian(numericCovfunHyp, x = 0, idx = idx)
    dim(dCov) <- c(nrow(sysDt), nrow(sysDt))
    numericGrad[idx] <- sum(diag(t(L1) %*% dCov %*% L1)) - sum(diag(t(L2) %*% dCov %*% L2))
  }
  expect_equal(analyticGrad[idcs], numericGrad[idcs])
})


test_that("grad_xCovX_wrtUnc produces correct result", {

  d <- runif(nrow(sysDt))
  analyticGrad <- as.vector(sysCompHandler$grad_xtCovx_wrtUnc(d, sysDt, ret.mat = TRUE))
  numericGrad <- rep(0, nrow(sysDt))
  idcs <- c(1,100)
  for (idx in idcs)
    numericGrad[idx] <- as.vector(jacobian(numericGradfunUnc, x = 0, idx = idx, z = d))
  expect_equal(analyticGrad[idcs], numericGrad[idcs])
  # subsetting
  analyticGrad <- as.vector(sysCompHandler$grad_xtCovx_wrtUnc(d, sysDt, idx = idcs, ret.mat = TRUE))
  expect_equal(analyticGrad, numericGrad)
})


test_that("gradLogLike wrt systematic uncertainties produces correct result", {

  setkey(sysDt, IDX)
  S <- sysCompHandler$map(exampleExpDt, sysDt, ret.mat = TRUE)
  D <- Diagonal(x = pmax(0.1 * exampleExpDt$DATA, 1)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)
  d <- exampleExpDt$DATA * runif(nrow(exampleExpDt), -0.1, 0.1)

  gradList <- gradLogLike(d, D, S, sysDt = sysDt, gpDt = gpDt, sysCompHandler = sysCompHandler)
  gradList$gradLogLike_wrtSysUnc[sysDt[!is.na(GPTYPE), IDX]] <- 0
  selIdcs <- order(abs(as.vector(gradList$gradLogLike_wrtSysUnc)), decreasing = TRUE)[sample(20, 5, replace = FALSE)]
  for (curIdx in selIdcs) {
    numRes <- as.vector(jacobian(numericGradfunLoglike, x = 0, idx = curIdx, d = d, D = D, S = S, parType = "sysUnc"))
    expect_equal(gradList$gradLogLike_wrtSysUnc[curIdx], numRes, tolerance = 5e-6)
  }

  # if a block in sysDt is associated with a GP, derivatives wrt UNC must be zero
  expect_true(all(as.vector(gradList$gradLogLike_wrtSysUnc) == 0 | is.na(sysDt$GPTYPE)))
})


test_that("gradLogLike wrt hyperparameters produces correct result", {

  S <- sysCompHandler$map(exampleExpDt, sysDt, ret.mat = TRUE)
  D <- Diagonal(x = pmax(0.1 * exampleExpDt$DATA, 1)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)
  d <- exampleExpDt$DATA * runif(nrow(exampleExpDt), -0.1, 0.1)

  system.time(gradList <- gradLogLike(d, D, S, sysDt = sysDt, gpDt = gpDt,
                                                statUncIdx = NULL, gpIdx = NULL, sysCompHandler = sysCompHandler))
  selIdcs <- order(abs(as.vector(gradList$gradLogLike_wrtHyp)), decreasing = TRUE)[1:10]
  for (curIdx in selIdcs) {
    numRes <- as.vector(jacobian(numericGradfunLoglike, x = 0, idx = curIdx, d = d, D = D, S = S, parType = "gpPar"))
    expect_equal(gradList$gradLogLike_wrtHyp[curIdx], numRes, tolerance = 5e-6)
  }
})


test_that("gradLogLike wrt statistical uncertainties produces correct result", {

  S <- sysCompHandler$map(exampleExpDt, sysDt, ret.mat = TRUE)
  D <- Diagonal(x = pmax(0.1 * exampleExpDt$DATA, 1)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)
  d <- exampleExpDt$DATA * runif(nrow(exampleExpDt), -0.1, 0.1)

  system.time(gradList <- gradLogLike(d, D, S, sysDt = sysDt, gpDt = gpDt,
                                      statUncIdx = NULL, gpIdx = NULL, sysCompHandler = sysCompHandler))
  selIdcs <- order(abs(as.vector(gradList$gradLogLike_wrtStatUnc)), decreasing = TRUE)[1:10]
  for (curIdx in selIdcs) {
    numRes <- as.vector(jacobian(numericGradfunLoglike, x = 0, idx = curIdx, d = d, D = D, S = S, parType = "statUnc"))
    expect_equal(gradList$gradLogLike_wrtStatUnc[curIdx], numRes, tolerance = 5e-6)
  }
})





