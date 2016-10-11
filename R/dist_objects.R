#' General components of Exponential Family distribution
#'
#' @param h (function) the measure
#' @param eta.from.theta (function) canonical parameter from mean parameter
#' @param theta.from.eta (function) mean parameter from canonical parameter
#' @param B.from.theta (function) value of log-normalizer for mean parameter
#' @param A.from.eta (function) value of log-normalizer for natural parameter in canonical parametrization
#' @param S (function) vector of sufficient statistics
ExpFam_dist <- function(
  h,
  eta.from.theta,
  theta.from.eta,
  B.from.theta,
  A.from.eta,
  S
){
  dst <- list(
    h = h,
    eta.from.theta = eta.from.theta,
    theta.from.eta = theta.from.eta,
    B.from.theta = B.from.theta,
    A.from.eta = A.from.eta,
    S = S
  )
  class(dst) <- "ExpFam_dist"
  dst
}

#' The generic density function for canonical parametrization
ExpFam_density_theta <-
  function(x, ...) UseMethod("ExpFam_density_theta")

#' The density function for canonical parametrization
ExpFam_density_theta.ExpFam_dist <- function(org.dist, theta){
  theta.bounded <- theta
  eta <- org.dist$eta.from.theta(theta.bounded)
  B.theta.bounded <- org.dist$B.from.theta(theta)
  dens <- function(x){
    org.dist$h(x) *
      exp( sum( eta * org.dist$S(x) ) - B.theta.bounded)
  }
  dens
}

#' The generic density function for mean parametrization
ExpFam_density_eta <- function(x, ...) UseMethod("ExpFam_density_eta")


#' The density function for mean parametrization
ExpFam_density_eta.ExpFam_dist <- function(org.dist, eta){
  eta.bounded <- eta
  A.eta.bounded <- org.dist$A.from.eta(eta.bounded)
  dens <- function(x){
    org.dist$h(x)*exp( sum( eta.bounded* org.dist$S(x)) - A.eta.bounded)
  }
  dens
}

