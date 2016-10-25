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
#'
#' \eqn{d(x|\theta) = h(x) \exp{(\eta(\theta)^T S(x) - B(\theta))}}
#'
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
    # ncol(eta) == org.dist$eta.dim
    dens <- function(x){
      org.dist$h(x) *
        exp( eta %*% org.dist$S(x) - B.theta.bounded)
    }
  }
  dens
}

#' The generic density function for mean parametrization
#' @export
ExpFam_density_eta <- function(x, ...) UseMethod("ExpFam_density_eta")


#' The density function for mean parametrization
#'
#' \eqn{d(x|\eta) = h(x) \exp{(\eta^T S(x) - A(\eta))}}
#'
#' @export
ExpFam_density_eta.ExpFam_dist <- function(org.dist, eta, num.opt = TRUE, verify = TRUE){
  ExpFam_density.ExpFam_dist_canonical(
    org.dist = org.dist,
    eta = eta,
    num.opt = num.opt,
    verify = verify
  )
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
  #moments of suff. stats
  if (!is.null(z$dA.deta)) {
    z$mean.S.eta <- function(eta) z$dA.deta(eta)
  }
  if (!is.null(z$d2A.deta2)) {
    z$var.S.eta <- function(eta) z$d2A.deta2(eta)
  }
  # loglik
  if (is.null(z$N.logLik.eta)) {
    z$N.logLik.eta <- function(x, eta) z$N.density.eta(x = x, eta = eta, log = TRUE)
  }
  if (is.null(z$N.logLik.theta)) {
    z$N.logLik.theta <- function(x, theta) z$N.density.theta(x = x, theta = theta, log = TRUE)
  }

  z
}



conjugate_distribution_eta.ExpFam_dist <- function(org.dist, h.c, A.from.eta.c) {
  #using natural parameters
  #in conjugate prior x <~~ \eta
  ExpFam_dist(
    h = h.c,
    eta.from.theta = NULL,
    theta.from.eta = NULL,
    B.from.theta = NULL,
    A.from.eta = A.from.eta.c,
    S = function(x) cbind(x, -org.dist$A.from.eta(x)),
    eta.dim = org.dist$eta.dim + 1
  )
}

tmptmp2 <- function() {
  bern.obj <- dist_Bernoulli()
  dd <- function(a, b, x) a*x - b*log(1 + exp(x))
  exp(x)^a/(1+exp(x))^b

  dd
}

tmptmp <- function() {
  bern.obj <- dist_Bernoulli()
  bern.obj$B.from.theta
  bern.obj$theta.from.eta
  bern.obj$stats.glm$mu.eta(2)   # dmu/deta
  bern.obj$A.from.eta
  str(bern.obj)
  bern.obj$eta.from.theta
  bern.obj$eta.from.theta(0.001)
  bern.obj$eta.from.theta(1-0.001)
  log(theta) - log(1 - theta)
  "1/theta + 1/(1-theta) = 1/((1-theta)*theta)"
  exp(a*(log(theta) - log(1 - theta)) -log(1 - theta)*b)

  beta.obj <- dist_Beta()

  bern.obj$d2A.deta2(2)
  bern.obj$d2A.deta2(bern.obj$eta.from.theta(0.4))
  bern.obj$stats.glm$mu.eta(bern.obj$N.eta.from.theta(1 - 0.01))
  dbeta(0.3, 2, 3)
  gamma(2) + gamma(3) - gamma(5)
  t2 <- log(gamma(2)) + log(gamma(3)) - log(gamma(5))
  t1 <- 2*log(0.3) + 3*log(0.7) # \eta^T S
  exp(t1 -t2)/(0.3*0.7)

  binomial()
  bern.obj$stats.glm$linkfun(0.5)
}
