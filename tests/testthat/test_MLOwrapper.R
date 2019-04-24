library(testthat)
context("MLO wrapper")

library(data.table)
library(jsonExforUtils)

# load example dataset

expDt <- copy(exampleExpDt)
subents <- copy(exampleSubents)

# create some systematic components

expHandler <- createSysCompHandler()

reacHandler <- createSysCompReacHandler(subents)
normHandler <- createSysCompNormHandler("REFDATA")
gpHandler <- createSysCompGPHandler()

# initialize data for handlers

reacHandler$addMap("pw", pwMap)
curEnergyGrid <- seq(0,31, length=30)

curReac <- "(26-FE-56(N,P)25-MN-56,,SIG)"
curEnGrid <- seq(0, 30, by = 1)  # seq(getThresEn(curReac, modDt), 32, by = 1)
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = c(1e4, 1e4, rep(20, length(curEnGrid)-2)),
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

normHandler$addSysUnc("EXPID", unique(expDt$EXPID), 0, 0.1)
gpHandler$addGP("REACEXP-N,P-01", 10, 5, 1e-4)


expHandler$addHandler(reacHandler)
expHandler$addHandler(normHandler)
expHandler$addGPHandler(gpHandler)

sysDt <- expHandler$createSysDt()
gpDt <- gpHandler$createGPDt()
gpDt[, IDX := seq_len(.N)]
setkey(sysDt, IDX, DIDX)
sysDt[EXPID=="REACEXP-N,P-01", EN := seq(1,30,length=.N)]


# necessary additional information

test_that("logLike correctly calculated for varied statistical and systematic components (statUnc/sysUnc)", {

  expDt[, REFDATA := runif(.N)]
  sysDt[, REFDATA := runif(.N)]

  curAdjStat <- sample(c(FALSE,TRUE), nrow(expDt), replace = TRUE)
  curAdjSys <- sample(c(FALSE,TRUE), nrow(sysDt), replace = TRUE)

  setkey(expDt, IDX)
  expDt[, ADJUSTABLE := curAdjStat]
  setkey(sysDt, IDX)
  sysDt[, ADJUSTABLE := curAdjSys]

  optfuns <- createMLOptimFuns()
  optfuns$setDts(expDt, sysDt, NULL, expHandler)

  setkey(expDt, IDX)
  setkey(sysDt, IDX)
  testx <- c(expDt[ADJUSTABLE == TRUE, DATA * runif(.N, 0.1, 0.2)],
             sysDt[ADJUSTABLE == TRUE, runif(.N, 0.1, 0.2)])

  D <- Diagonal(x = getDt_UNC(expDt)^2)
  diag(D)[curAdjStat] <- testx[seq_len(sum(expDt$ADJUSTABLE))]^2
  S <- expHandler$map(expDt, sysDt, ret.mat = TRUE)
  P <- expHandler$cov(sysDt, ret.mat = TRUE)
  diag(P)[curAdjSys] <- testx[sum(expDt$ADJUSTABLE) + seq_len(sum(sysDt$ADJUSTABLE))]^2
  d <- getDt_DATA(expDt) - getDt_REFDATA(expDt) -
    S %*% (getDt_DATA(sysDt) - getDt_REFDATA(sysDt))


  res1 <- logLike(d, D, S, P)
  res2 <- optfuns$logLike(testx)
  expect_equal(res1, res2)
})


test_that("modified UNC element in sysDt masked by a GP should not affect result", {

  gpSysDt <- copy(sysDt)
  gpHandler$updateSysDt(gpSysDt)
  gpDt[, ADJUSTABLE := FALSE]

  expDt[, REFDATA := runif(.N)]
  sysDt[, REFDATA := runif(.N)]

  curAdjStat <- sample(c(FALSE,TRUE), nrow(expDt), replace = TRUE)
  curAdjSys <- sample(c(FALSE,TRUE), nrow(gpSysDt), replace = TRUE)

  setkey(expDt, IDX)
  expDt[, ADJUSTABLE := curAdjStat]
  setkey(gpSysDt, IDX)
  gpSysDt[, ADJUSTABLE := curAdjSys]
  isGP <- !is.na(gpSysDt$GPTYPE)

  res1 <- expHandler$cov(gpSysDt, gpDt, ret.mat = TRUE)
  gpSysDt[!is.na(GPTYPE), UNC := runif(.N)]
  res2 <- expHandler$cov(gpSysDt, gpDt, ret.mat = TRUE)
  expect_equal(res1, res2)
})


test_that("logLike correctly calculated if GP component included", {

  gpSysDt <- copy(sysDt)
  gpHandler$updateSysDt(gpSysDt)
  gpDt <- copy(gpDt)

  expDt[, REFDATA := runif(.N)]
  sysDt[, REFDATA := runif(.N)]

  curAdjStat <- sample(c(FALSE,TRUE), nrow(expDt), replace = TRUE)
  curAdjSys <- sample(c(FALSE,TRUE), nrow(gpSysDt), replace = TRUE)

  setkey(expDt, IDX)
  expDt[, ADJUSTABLE := curAdjStat]
  setkey(gpSysDt, IDX)
  gpSysDt[, ADJUSTABLE := curAdjSys]
  isGP <- !is.na(gpSysDt$GPTYPE)
  setkey(gpDt, IDX)
  gpDt[, ADJUSTABLE := FALSE]

  optfuns <- createMLOptimFuns()
  optfuns$setDts(expDt, gpSysDt, gpDt, expHandler)

  setkey(expDt, IDX)
  setkey(gpSysDt, IDX)
  testx <- c(expDt[ADJUSTABLE == TRUE, DATA * runif(.N, 0.1, 0.2)],
             gpSysDt[ADJUSTABLE == TRUE, runif(.N, 0.1, 0.2)])

  D <- Diagonal(x = getDt_UNC(expDt)^2)
  diag(D)[curAdjStat] <- testx[seq_len(sum(expDt$ADJUSTABLE))]^2
  S <- expHandler$map(expDt, gpSysDt, ret.mat = TRUE)
  P <- expHandler$cov(gpSysDt, gpDt, ret.mat = TRUE)

  setkey(gpSysDt, IDX)
  diag(P)[curAdjSys & !isGP] <- testx[sum(expDt$ADJUSTABLE) + which((gpSysDt$ADJUSTABLE & !isGP)[gpSysDt$ADJUSTABLE])]^2
  d <- getDt_DATA(expDt) - getDt_REFDATA(expDt) -
    S %*% (getDt_DATA(gpSysDt) - getDt_REFDATA(gpSysDt))


  res1 <- logLike(d, D, S, P)
  res2 <- optfuns$logLike(testx)
  expect_equal(res1, res2)
})



test_that("logLike correctly calculated with varied GP component", {

  expDt <- copy(expDt)
  gpSysDt <- copy(sysDt)
  gpHandler$updateSysDt(gpSysDt)
  gpDt <- copy(gpDt)

  expDt[, REFDATA := runif(.N)]
  sysDt[, REFDATA := runif(.N)]

  curAdjStat <- sample(c(FALSE,TRUE), nrow(expDt), replace = TRUE)
  curAdjSys <- sample(c(FALSE,TRUE), nrow(gpSysDt), replace = TRUE)
  curAdjGp <- (gpDt$PARNAME %in% c("sigma","len"))[order(gpDt$IDX)]

  setkey(expDt, IDX)
  expDt[, ADJUSTABLE := curAdjStat]
  setkey(gpSysDt, IDX)
  gpSysDt[, ADJUSTABLE := curAdjSys]
  isGP <- !is.na(gpSysDt$GPTYPE)
  setkey(gpDt, IDX)
  gpDt[, ADJUSTABLE := curAdjGp]

  optfuns <- createMLOptimFuns()
  optfuns$setDts(expDt, gpSysDt, gpDt, expHandler)

  setkey(expDt, IDX)
  setkey(gpSysDt, IDX)

  expDt[ADJUSTABLE == TRUE, UNC := DATA * runif(.N, 0.1, 0.2)]
  gpSysDt[ADJUSTABLE == TRUE, UNC := runif(.N, 0.1, 0.2)]
  gpDt[ADJUSTABLE == TRUE, PARVAL := runif(.N, 0.1, 0.2)]

  testx <- c(expDt[ADJUSTABLE == TRUE, UNC],
             gpSysDt[ADJUSTABLE == TRUE, UNC],
             gpDt[ADJUSTABLE == TRUE, PARVAL])

  D <- Diagonal(x = getDt_UNC(expDt)^2)
  diag(D)[curAdjStat] <- testx[seq_len(sum(expDt$ADJUSTABLE))]^2
  S <- expHandler$map(expDt, gpSysDt, ret.mat = TRUE)
  P <- expHandler$cov(gpSysDt, gpDt, ret.mat = TRUE)

  d <- getDt_DATA(expDt) - getDt_REFDATA(expDt) -
    S %*% (getDt_DATA(gpSysDt) - getDt_REFDATA(gpSysDt))

  res1 <- logLike(d, D, S, P)
  res2 <- optfuns$logLike(testx)
  expect_equal(res1, res2)
})



test_that("gradLogLike correctly calculated with varied statUnc, sysUnc, and gpHyps", {

  expDt <- copy(expDt)
  gpSysDt <- copy(sysDt)
  gpHandler$updateSysDt(gpSysDt)
  gpDt <- copy(gpDt)

  expDt[, REFDATA := runif(.N)]
  gpSysDt[, REFDATA := runif(.N)]

  curAdjStat <- sample(c(FALSE,TRUE), nrow(expDt), replace = TRUE)
  curAdjSys <- sample(c(FALSE,TRUE), nrow(gpSysDt), replace = TRUE)
  curAdjGp <- (gpDt$PARNAME %in% c("sigma","len"))[order(gpDt$IDX)]
  statUncIdx <- which(curAdjStat)
  sysUncIdx <- which(curAdjSys)
  gpIdx <- which(curAdjGp)

  setkey(expDt, IDX)
  expDt[, ADJUSTABLE := curAdjStat]
  setkey(gpSysDt, IDX)
  gpSysDt[, ADJUSTABLE := curAdjSys]
  isGP <- !is.na(gpSysDt$GPTYPE)
  setkey(gpDt, IDX)
  gpDt[, ADJUSTABLE := curAdjGp]

  optfuns <- createMLOptimFuns()
  optfuns$setDts(expDt, gpSysDt, gpDt, expHandler)

  setkey(expDt, IDX)
  setkey(gpSysDt, IDX)

  expDt[ADJUSTABLE == TRUE, UNC := DATA * runif(.N, 0.1, 0.2)]
  gpSysDt[ADJUSTABLE == TRUE, UNC := runif(.N, 0.1, 0.2)]
  gpDt[ADJUSTABLE == TRUE, PARVAL := runif(.N, 0.1, 0.2)]

  testx <- c(expDt[ADJUSTABLE == TRUE, UNC],
             gpSysDt[ADJUSTABLE == TRUE, UNC],
             gpDt[ADJUSTABLE == TRUE, PARVAL])

  D <- Diagonal(x = getDt_UNC(expDt)^2)
  diag(D)[curAdjStat] <- testx[seq_len(sum(expDt$ADJUSTABLE))]^2
  S <- expHandler$map(expDt, gpSysDt, ret.mat = TRUE)
  P <- expHandler$cov(gpSysDt, gpDt, ret.mat = TRUE)

  d <- getDt_DATA(expDt) - getDt_REFDATA(expDt) -
    S %*% (getDt_DATA(gpSysDt) - getDt_REFDATA(gpSysDt))

  rawRes1 <- gradLogLike(d, D, S,
                      statUncIdx = statUncIdx,
                      sysUncIdx = sysUncIdx,
                      gpIdx = gpIdx, sysCompHandler = expHandler,
                      sysDt = gpSysDt, gpDt = gpDt)
  res1 <- with(rawRes1, c(gradLogLike_wrtStatUnc[statUncIdx],
                          gradLogLike_wrtSysUnc[sysUncIdx],
                          gradLogLike_wrtHyp[gpIdx]))

  res2 <- optfuns$gradLogLike(testx)
  expect_equal(res1, res2)
})












