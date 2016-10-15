#' General components of Exponential Family distribution
#'
#' @param h (function) the measure
#' @param eta.from.theta (function) canonical parameter from mean parameter
#' @param theta.from.eta (function) mean parameter from canonical parameter
#' @param B.from.theta (function) value of log-normalizer for mean parameter
#' @param A.from.eta (function) value of log-normalizer for natural parameter in canonical parametrization
#' @param S (function) vector of sufficient statistics
#' @param eta.dim (numeric), number of parameters of the exponential family being defined
ExpFam_dist <- function(
  h,
  eta.from.theta,
  theta.from.eta,
  B.from.theta,
  A.from.eta,
  S,
  eta.dim
){
  dst <- list(
    h = h,
    eta.from.theta = eta.from.theta,
    theta.from.eta = theta.from.eta,
    B.from.theta = B.from.theta,
    A.from.eta = A.from.eta,
    S = S,
    eta.dim = eta.dim
  )
  class(dst) <- "ExpFam_dist"
  dst
}

#' The generic density function for canonical parametrization
#' @export
ExpFam_density_theta <-
  function(x, ...) UseMethod("ExpFam_density_theta")

#' The density function for canonical parametrization
#' @export
ExpFam_density_theta.ExpFam_dist <- function(org.dist, theta){
  theta.bounded <- theta
  eta <- org.dist$eta.from.theta(theta.bounded)
  B.theta.bounded <- org.dist$B.from.theta(theta)
  if (org.dist$eta.dim == 1) {
    dens <- function(x){
      org.dist$h(x) *
        exp( eta * org.dist$S(x) - B.theta.bounded)
    }
  } else {
    dens <- function(x){
      org.dist$h(x) *
        exp( sum( eta * org.dist$S(x) ) - B.theta.bounded)
    }
    warning("Buggy code: improper summation for eta.dim > 1")
  }
  dens
}

#' The generic density function for mean parametrization
#' @export
ExpFam_density_eta <- function(x, ...) UseMethod("ExpFam_density_eta")


#' The density function for mean parametrization
#' @export
ExpFam_density_eta.ExpFam_dist <- function(org.dist, eta){
  eta.bounded <- eta
  A.eta.bounded <- org.dist$A.from.eta(eta.bounded)
  if (org.dist$eta.dim == 1) {
    dens <- function(x){
      org.dist$h(x)*exp( eta.bounded* org.dist$S(x) - A.eta.bounded)
    }
  } else {
    dens <- function(x){
      org.dist$h(x)*exp( sum( eta.bounded* org.dist$S(x)) - A.eta.bounded)
    }
    warning("Buggy code: improper summation for eta.dim > 1")
  }
  dens
}

#' The generic property of having a GLM-compatible object
#' @export
has_stat_GLM <- function(x, ...) UseMethod("has_stat_GLM")

#' @export
has_stat_GLM.ExpFam_dist <- function(org.dist) {
  !is.null(org.dist$stats.glm)
}

#' @export
N_autogenerate_from_stat_GLM <- function(x, ...) UseMethod("N_autogenerate_from_stat_GLM")


#' @export
N_autogenerate_from_stat_GLM.ExpFam_dist <- function(z) {
  if (has_stat_GLM(z)) {
    z$N.eta.from.theta <- function(theta) z$stats.glm$linkfun(theta)
    z$N.theta.from.eta <- function(eta) z$stats.glm$linkinv(eta)
    z$N.variance.from.mean <- function(mu) z$stats.glm$variance(mu)
  }
  z
}

#' @export
N_autogenerate_functions <- function(x, ...) UseMethod("N_autogenerate_functions")


#' @export
N_autogenerate_functions.ExpFam_dist <- function(z) {
  if (is.null(z$N.density.eta) && !is.null(z$N.density.theta)) {
    #autodef N.density.eta
    z$N.density.eta <- function(x, eta, log = FALSE) z$N.density.theta(x, theta = z$theta.from.eta(eta), log = log)
  }
  if (is.null(z$N.density.theta) && !is.null(z$N.density.eta)) {
    #autodef N.density.eta
    z$N.density.theta <- function(x, theta, log = FALSE) z$N.density.eta(x, eta = z$eta.from.theta(theta), log = log)
  }
  if (is.null(z$N.logLik.eta)) {
    z$N.logLik.eta <- function(x, eta) z$N.density.eta(x = x, eta = eta, log = TRUE)
  }
  if (is.null(z$N.logLik.theta)) {
    z$N.logLik.theta <- function(x, theta) z$N.density.theta(x = x, theta = theta, log = TRUE)
  }
  z
}
