#' Components of Exponential Family distribution
#'
#' @param h (function) the measure
#' @param A.from.eta (function) value of log-normalizer for natural parameter in canonical parametrization
#' @param S (function) vector of sufficient statistics
#' @param eta.dim (numeric), number of parameters of the exponential family being defined
#' @param A.grad (function) gradient of the function A.from.eta
#' @param A.hess (function) hessian of the function A.from.eta
#' @param eta.in.domain
#' @export
ExpFam_dist_canonical <- function(
  h,
  A.from.eta,
  S,
  eta.dim,
  A.grad,
  A.hess,
  eta.in.domain
){
  dst <- list(
    h = h,
    A.from.eta = A.from.eta,
    S = S,
    eta.dim = eta.dim,
    A.grad = A.grad,
    A.hess = A.hess,
    eta.in.domain = eta.in.domain
  )
  class(dst) <- "ExpFam_dist_canonical"
  dst
}

#' @export
ExpFam_param <- function(
  org.dist,
  eta.from.theta,
  theta.from.eta,
  B.from.theta,
  B.grad,
  B.hess,
  theta.in.domain,
  param.type
){
  dst <- list(
    org.dist = org.dist,
    eta.from.theta = eta.from.theta,
    theta.from.eta = theta.from.eta,
    B.from.theta = B.from.theta,
    B.grad = B.grad,
    B.hess = B.hess,
    theta.in.domain = theta.in.domain,
    param.type = param.type
  )
  class(dst) <- "ExpFam_param"
  dst
}

#' @export
ExpFam_dist_ext <- function(
  canonical.param,
  parametrisations
) {

  z <- parametrisations
  names(z) <- lapply(z, function(x) x$param.type) %>% unlist()
  z <- c(z, list(canonical = canonical.param))
  class(z) <- "ExpFam_dist_ext"
  z
}




#' The generic density function for canonical parametrization
#' @export
ExpFam_density <- function(x, ...) UseMethod("ExpFam_density")


#' The density function for canonical parametrization
#'
#' \eqn{d(x|\eta) = h(x) \exp{(\eta^T S(x) - A(\eta))}}
#'
#' @export
ExpFam_density.ExpFam_dist_canonical <- function(org.dist, eta, num.opt = TRUE, verify = TRUE){
  if (verify && !functor_eval_1(org.dist$eta.in.domain, default.val = TRUE, eta = eta)) {
    stop("Eta not in domain")
  }
  if (num.opt && !is.null(org.dist$N.density.eta)) {
    eta.bound <- eta
    dens <- function(x) org.dist$N.density.eta(x = x, eta = eta.bound, log = FALSE)
  } else {
    if (org.dist$eta.dim == 1) {
      eta.bound <- eta
      A.eta.bound <- f_A_from_eta(org.dist)(eta.bound)
      dens <- function(x){
        org.dist$h(x)*exp( eta.bound* org.dist$S(x) - A.eta.bound)
      }
    } else {
      # ncol(eta) == org.dist$eta.dim
      eta.bound <- matrix(eta, nrow = 1)
      A.eta.bound <- org.dist$A.from.eta(eta.bound)
      dens <- function(x){
        org.dist$h(x)*exp( org.dist$S(x) %*% t(eta.bound)  - A.eta.bound)[,1]
      }
    }
  }
  dens
}

#' @export
ExpFam_density.ExpFam_param <- function(par.dist, theta, num.opt = TRUE, verify = TRUE) {
  if (verify && !functor_eval_1(par.dist$theta.in.domain, default.val = TRUE, theta = theta)) {
    stop("Theta not in domain")
  }
  theta.bound <- theta
  eta <- f_eta_from_theta(par.dist, num.opt = num.opt)(theta.bound)
  B.theta.bound <- functor_eval_first_not_null(
    c(
      f_B_from_theta(par.dist, num.opt = num.opt, error.cond = FALSE),
      function(theta){
        f_A_from_eta(par.dist$org.dist)(eta)
      }
    ),
    theta = theta
  )
  if (par.dist$org.dist$eta.dim == 1) {
    dens <- function(x){
      par.dist$org.dist$h(x) *
        exp( eta * par.dist$org.dist$S(x) - B.theta.bound)
    }
  } else {
    # ncol(eta) == org.dist$eta.dim
    dens <- function(x){
      par.dist$org.dist$h(x) *
        exp( eta %*% par.dist$org.dist$S(x) - B.theta.bound)
    }
  }
  dens

}

#' @export
has_stat_GLM.ExpFam_param <- function(param.dist) {
  !is.null(param.dist$stats.glm)
}


#' @export
N_autogenerate_from_stat_GLM.ExpFam_param <- function(z) {
  if (has_stat_GLM(z)) {
    z$N.eta.from.theta <- function(theta) z$stats.glm$linkfun(theta)
    z$N.theta.from.eta <- function(eta) z$stats.glm$linkinv(eta)
    if (z$param.type == "mean") {
      z$N.theta.from.eta.grad <- function(eta) z$stats.glm$mu.eta(eta)
      z$N.variance.from.theta <- function(theta) z$stats.glm$variance(theta)
    }
  }
  z
}

#' @export
conjugate_distribution_proto.ExpFam_dist_canonical <- function(org.dist) {
  #using natural parameters
  #in conjugate prior x <~~ \eta
  ExpFam_dist_canonical(
    h = org.dist$conjugate.prior.elems_raw$h,
    A.from.eta = org.dist$conjugate.prior.elems_raw$A.from.eta,
    S = function(x) cbind(x, -org.dist$A.from.eta(x)),
    eta.dim = org.dist$eta.dim + 1,
    A.grad = org.dist$conjugate.prior.elems_raw$A.grad,
    A.hess = org.dist$conjugate.prior.elems_raw$A.hess,
    eta.in.domain = org.dist$conjugate.prior.elems_raw$eta.in.domain
  )
}

tmptmptmp <- function(){
  tmp <- ExpFam_dist_Bernoulli2()
  ccc <- conjugate_distribution_proto.ExpFam_dist_canonical(tmp$par.canonical)
  ExpFam_density(ccc, c(1,2), verify = FALSE)(-1)
  s.tmp <- seq(from = -10, to = 10, by = 1)
  ExpFam_density(ccc, c(1,2), verify = FALSE)(s.tmp)
  dbeta(tmp$par.mean$theta.from.eta(s.tmp), 1.01, 1.01)
}
