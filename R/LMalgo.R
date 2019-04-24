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


#' Levenberg-Marquardt algorithm with prior
#'
#' @param fn model function
#' @param jac Jacobian of model function
#' @param pinit initial model parameters
#' @param p0 prior mean of model parameters
#' @param P0 prior covariance matrix of model parameters
#' @param yexp experimental data points
#' @param D only diagonal contribution to experimental covariance matrix
#' @param S mapping from systematic components to observations
#' @param X prior covariance matrix for systematic components
#' @param control control arguments, see \code{details}
#'
#' @details
#' The control argument is a list that can contain
#' \itemize{
#'   \item maximum number of iterations \code{maxit}
#'   \item absolute tolerance required for convergence \code{abstol}
#'   \item relative tolerance required for convergence \code{reltol}
#' }
#'
#' @return
#' list with elements \code{par}, \code{fn}, \code{jac}, \code{counts}, \code{value}
#' @export
#'
LMalgo <- function(fn, jac, pinit, p0, P0, yexp, D, S, X,
                   lower = -Inf, upper = Inf, logger = NULL,
                   control = list(maxit = 10, abstol = 1e-3, reltol = 1e-3, mu = NULL)) {

  parDefaults <- list(
    control = list(maxit = 10, abstol = 1e-3, reltol = 1e-3)
  )
  control <- modifyList(parDefaults$control, control)

  isDebugOn <- TRUE
  debuginfo <- function(str) {
    if (isDebugOn)
      cat(str, "\n")
  }

  # initialization
  tau <- 1e-3  # factor to determine damping factor
  DD <- Diagonal(x = rep(1, length(pinit)))
  firstLoop <- TRUE

  # speeds up inversions of the experimental covariance matrix
  cholZ <- makeCholZ(D, S, X)
  breakCounter <- 0

  invP0 <- solve(P0)
  pref <- as.vector(pinit)
  J <- jac(as.vector(pref))
  fref <- fn(as.vector(pref))

  logBuffer <- list()
  i <- 0  # counter
  while (i < control$maxit && breakCounter < 3) {
    i <- i + 1

    print(breakCounter)
    # propose new parameter set and calculate function value
    tJinvBJ <- forceSymmetric(mult_xt_invCov_x(J, D, S, X, cholZ = cholZ))
    invP1 <- invP0 + tJinvBJ

    if (firstLoop) {
      mu <- if (is.null(control$mu)) tau * max(diag(invP1)) else control$mu
      firstLoop <- FALSE
    }

    m <- t(J) %*% mult_invCov_x(yexp - fref, D, S, X, cholZ = cholZ) +
      invP0 %*% (p0 - pref)
    pprop <- as.vector(pref + solve(invP1 + mu * DD, m))

    # project on box constraints
    pprop_unconstr <- pprop
    pprop[pprop < lower] <- lower[pprop < lower]
    pprop[pprop > upper] <- upper[pprop > upper]

    fprop <- fn(as.vector(pprop))
    fprop_approx <- fref + J %*% (pprop - pref)
    fprop_approx_unconstr <- fref + J %*% (pprop_unconstr - pref)

    # calculate true objective function and approximation thereof
    dpriorRef <- pref - p0
    dpriorProp <- pprop - p0
    dpriorProp_unconstr <- pprop_unconstr - p0

    LpriorRef <- as.vector(crossprod(dpriorRef, invP0 %*% dpriorRef))
    LpriorProp <- as.vector(crossprod(dpriorProp, invP0 %*% dpriorProp))
    LpriorProp_unconstr <- as.vector(crossprod(dpriorProp_unconstr, invP0 %*% dpriorProp_unconstr))

    Lref <- as.vector(mult_xt_invCov_x(yexp - fref, D, S, X, cholZ = cholZ)) + LpriorRef
    Lprop <- as.vector(mult_xt_invCov_x(yexp - fprop, D, S, X, cholZ = cholZ)) + LpriorProp
    Lprop_approx <- as.vector(mult_xt_invCov_x(yexp - fprop_approx, D, S, X, cholZ = cholZ)) + LpriorProp
    Lprop_approx_unconstr <- as.vector(mult_xt_invCov_x(yexp - fprop_approx_unconstr, D, S, X, cholZ = cholZ)) + LpriorProp_unconstr

    stopifnot(Lprop_approx_unconstr < Lref)

    # adjust step size
    # Marquardt 1963
    gain <- ((Lref - Lprop)+1e-10) / (abs(Lref - Lprop_approx)+1e-10)  # +1e-10 to avoid problems if Lprop_approx == Lref
                                                                       # abs in denominator important because of
                                                                       # projecting back into feasible region
                                                                       # and then Lprop_approx < Lref not guaranteed anymore
    if (gain < 0.25) mu <- mu * 2
    else if (gain > 0.75) mu <- mu / 3

    # check break conditions
    if (abs(Lprop - Lref) / abs(Lref) < control$reltol ||
        abs(Lprop - Lref) < control$abstol) {
      breakCounter <- breakCounter + 1
    } else {
      breakCounter <- 0
    }

    # print status information
    logBuffer <- list(iteration = i, mu = mu, gain = gain, pref = pref, fref = fref, Lref = Lref,
                      pprop = pprop, fprop = fprop, Lprop = Lprop, fprop_approx = fprop_approx, Lprop_approx = Lprop_approx,
                      Jref = J, p0 = p0, P0 = P0, yexp = yexp, D = D, S = S, X = X)
    if (is.function(logger)) logger(logBuffer)

    # accept if proposed parameter set better than old one
    if (Lprop < Lref) {
      pref <- pprop
      fref <- fprop
      J <- jac(as.vector(pprop))
      Lref <- Lprop
    }

  }
  list(par = pref,
       fn = fref,
       jac = J,
       counts = i,
       value = Lref)
}









