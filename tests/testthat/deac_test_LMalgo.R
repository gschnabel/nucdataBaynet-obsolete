context("LM algorithm")

expdata <- as.data.frame(
  rbind(c(2.513400000000E+00,  0.000000000000E+00),
        c(2.044333373291E+00,  5.000000000000E-02),
        c(1.668404436564E+00,  1.000000000000E-01),
        c(1.366418021208E+00,  1.500000000000E-01),
        c(1.123232487372E+00,  2.000000000000E-01),
        c(9.268897180037E-01,  2.500000000000E-01),
        c(7.679338563728E-01,  3.000000000000E-01),
        c(6.388775523106E-01,  3.500000000000E-01),
        c(5.337835317402E-01,  4.000000000000E-01),
        c(4.479363617347E-01,  4.500000000000E-01),
        c(3.775847884350E-01,  5.000000000000E-01),
        c(3.197393199326E-01,  5.500000000000E-01),
        c(2.720130773746E-01,  6.000000000000E-01),
        c(2.324965529032E-01,  6.500000000000E-01),
        c(1.996589546065E-01,  7.000000000000E-01),
        c(1.722704126914E-01,  7.500000000000E-01),
        c(1.493405660168E-01,  8.000000000000E-01),
        c(1.300700206922E-01,  8.500000000000E-01),
        c(1.138119324644E-01,  9.000000000000E-01),
        c(1.000415587559E-01,  9.500000000000E-01),
        c(8.833209084540E-02,  1.000000000000E+00),
        c(7.833544019350E-02,  1.050000000000E+00),
        c(6.976693743449E-02,  1.100000000000E+00),
        c(6.239312536719E-02,  1.150000000000E+00)))

fn <- function(beta) {
  E <- expdata[,2]
  beta[1]*exp(-beta[2]*E) + beta[3]*exp(-beta[4]*E) + beta[5]*exp(-beta[6]*E)
}

gr <- function(beta) {
  E <- expdata[,2]
  jacobian(fn, x = beta)

}

startpar <- c(0.5, 0.7, 3.6, 4.2, 4, 6.3)
truepar <- c(9.5100000027E-02, 1.0000000001E+00,
             8.6070000013E-01, 3.0000000002E+00,
             1.5575999998E+00, 5.0000000001E+00)

# simple case

dontrun <- function() {
  D <- Diagonal(x = rep(1, nrow(expdata)))
  S <- Diagonal(x = rep(1, nrow(expdata)))
  X <- Diagonal(x = rep(1, nrow(expdata)))

  p0 <- startpar
  P0 <- Diagonal(x = rep(1e4, length(p0)))
  yexp <- expdata[,1]

  pfit <- LMalgo(fn, gr, runif(6), p0, P0, yexp, D, S, X, control = list(maxiter=20))


  E <- 1:10
  fn <- function(x) {
    exp(as.vector(x)*E/10)
  }

  gr <- function(x) {
    as(cbind((E/10) * exp(x*E/10)), "sparseMatrix")
  }

  yexp <- fn(1.5)
  p0 <- 3
  P0 <- Diagonal(x = 100)
  D <- Diagonal(x = rep(1, length(yexp)))
  S <- Diagonal(x = rep(1, length(yexp)))
  X <- Diagonal(x = rep(1, length(yexp)))

  LMalgo(fn, gr, 10, p0, P0, yexp, D, S, X, lower = 1.7, control = list(maxit=20))

  plot(E, fn(5))
}




