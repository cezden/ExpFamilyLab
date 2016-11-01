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


#' The generic function for selecting a parametrization
#' @export
ExpFam_bind_parametrization <- function(x, ...) UseMethod("ExpFam_bind_parametrization")

#' ExpFam_dist with selected parametrisation
#' @param dist.obj object of the class \code{\link{ExpFam_dist_ext}}
#' @param param.name (\code{character}) the name of the parametrisation used as a theta
#' @export
ExpFam_bind_parametrization.ExpFam_dist_ext <- function(dist.obj, param.name, num.opt = TRUE){
  if (is.null(dist.obj[[param.name]])) {
    stop(paste0("Unknown parametrization '", param.name, "'"))
  }
  par.theta <- dist.obj[[param.name]]
  par.eta <- dist.obj[["canonical"]]
  z <- ExpFam_dist(
    h = par.eta$h,
    eta.from.theta = f_eta_from_theta(par.theta, num.opt = num.opt),
    theta.from.eta = f_theta_from_eta(par.theta, num.opt = num.opt),
    B.from.theta = f_B_from_theta(par.theta, num.opt = num.opt),
    A.from.eta = f_A_from_eta(par.eta, num.opt = num.opt),
    S = par.eta$S,
    eta.dim = par.eta$eta.dim
  )
  if (has_stat_GLM(par.theta)) {
    z$stats.glm <- par.theta$stats.glm
    z <- N_autogenerate_from_stat_GLM(z)
  }
  z
}

#' @export
ExpFam_reparametrize <- function(x, ...) UseMethod("ExpFam_reparametrize")

#' Exponential Family Reparametrisation
#' @param param.obj object of class \code{\link{ExpFam_param}}
#' @param reparam.obj reparametrisation descriptor,
#'  object of class \code{\link{ExpFam_dist_reparam}}
#' @param param.type the name for resulting parametrisation
#' @export
ExpFam_reparametrize.ExpFam_param <- function(param.obj, reparam.obj, param.type){

  eta.from.theta <- reparametrize_function_with_proto(
    f1 = param.obj$eta.from.theta,
    f2 = reparam.obj$y.from.x,
    res.proto = function(theta) NULL
  )
  theta.from.eta <- reparametrize_function_with_proto(
    f1 = reparam.obj$x.from.y,
    f2 = param.obj$theta.from.eta,
    res.proto = function(eta) NULL
  )
  B.from.theta <- reparametrize_function_with_proto(
    f1 = param.obj$B.from.theta,
    f2 = reparam.obj$y.from.x,
    res.proto = function(theta) NULL
  )
  B.grad <- reparametrize_function_grad_with_proto(
    f.from.z = param.obj$B.from.theta,
    f.from.z.grad = param.obj$B.grad,
    z.from.x = reparam.obj$y.from.x,
    z.from.x.grad = reparam.obj$y.grad,
    res.proto = function(theta) NULL
  )

  B.hess <- reparametrize_function_hess_with_proto(
    f.from.z = param.obj$B.from.theta,
    f.from.z.grad = param.obj$B.grad,
    f.from.z.hess = param.obj$B.hess,
    z.from.x = reparam.obj$y.from.x,
    z.from.x.grad = reparam.obj$y.grad,
    z.from.x.hess = reparam.obj$y.hess,
    res.proto = function(theta) NULL,
    simplify = TRUE
  )

  theta.in.domain <- function(theta)
    reparam.obj$x.in.domain(theta) &&
    param.obj$theta.in.domain(reparam.obj$y.from.x(theta))

  ExpFam_param(
    org.dist = param.obj$org.dist,
    eta.from.theta = eta.from.theta,
    theta.from.eta = theta.from.eta,
    B.from.theta = B.from.theta,
    B.grad = B.grad,
    B.hess = B.hess,
    theta.in.domain = theta.in.domain,
    param.type = param.type
  )
}
