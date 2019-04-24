context("marginal likelihood related functions")

library(testthat)
library(mvtnorm)
library(Matrix)

S <- as(matrix(runif(9*3), 9, 3), "sparseMatrix")
D <- Diagonal(x = runif(9,0.1,1)^2)
P <- matrix(runif(9), 3, 3)
P <- as(t(P) %*% P, "sparseMatrix")
x <- runif(9)


test_that("mult_invCov_x works correctly", {

  res1 <- mult_invCov_x(x, D, S, P)
  res2 <- solve(S %*% P %*% t(S) + D, x)
  expect_equal(res1, res2)
})


test_that("mult_xt_invCov_x works correctly", {

  res1 <- mult_xt_invCov_x(x, D, S, P)
  res2 <- x %*% solve(S %*% P %*% t(S) + D, x)
  expect_equal(res1, res2)
})


test_that("chisquare function works correctly", {

  expdata <- runif(nrow(D))
  moddata <- runif(nrow(D))
  d <- expdata - moddata
  res1 <- chisquare(d, D, S, P)
  res2 <- as.vector(t(d) %*% solve(S %*% P %*% t(S) + D) %*% d)
  expect_equal(res1, res2)
})


test_that("logDetCov function works correctly", {

  res1 <- logDetCov(D, S, P)
  res2 <- as.vector(determinant(S%*%P%*%t(S) + D)$modulus)
  expect_equal(res1, res2)
})


test_that("logLike function works correctly", {

  d <- runif(nrow(D))
  res1 <- logLike(d, D, S, P)
  res2 <- dmvnorm(d, rep(0, length(d)), as.matrix(D + S %*% P %*% t(S)), log = TRUE)
  expect_equal(res1, res2)
})








