#' Exponential Family Reparametrization
#'
#' The class describing 1-1 mapping of the parameters of the exponential family,
#' consisting of:
#' \itemize{
#'   \item{the mapping and inverse mapping,}
#'   \item{gradients,}
#'   \item{Hessians,}
#'   \item{domains characteristic functions}
#' }
#' All components are assumed to be vectorisable,
#' when dimension > 1 the mapped points are assumed to be passed as rows of the matrix.
#'
#' @param y.from.x (\code{function}) y = f(x)
#' @param x.from.y (\code{function}) x = f^{-1}(y) = finv(y)
#' @param y.grad (\code{function}) f'(x)
#' @param x.grad (\code{function}) finv'(y)
#' @param y.hess (\code{function}) f''(x)
#' @param x.hess (\code{function}) finv''(y)
#' @param x.in.domain (\code{function})
#' @param y.in.domain (\code{function}) returns \code{TRUE} iff all
#' @export
ExpFam_dist_reparam <- function(
  y.from.x,
  x.from.y,
  y.grad,
  x.grad,
  y.hess,
  x.hess,
  x.in.domain,
  y.in.domain
){
  z <- list(
    y.from.x = y.from.x,
    x.from.y = x.from.y,
    y.grad = y.grad,
    x.grad = x.grad,
    y.hess = y.hess,
    x.hess = x.hess,
    x.in.domain = x.in.domain,
    y.in.domain = y.in.domain

  )
  class(z) <- "ExpFam_dist_reparam"
  z
}

#' @export
N_autogenerate_from_stat_GLM.ExpFam_dist_reparam <- function(z){
  if (!is.null(z$stats.glm.link)) {
    z$N.y.from.x <- function(x) z$stats.glm.link$linkfun(x)
    z$N.x.from.y <- function(y) z$stats.glm.link$linkinv(y)
    z$N.x.grad <- function(y) z$stats.glm.link$mu.eta(y)
  }
  z
}

#' Inverse reparametrisation (1D)
#' @export
Reparam_Inverse <- function(){
  ExpFam_dist_reparam(
    y.from.x = function(x) 1/x,
    x.from.y = function(y) 1/y,
    y.grad = function(x) -1/x ^ 2,
    x.grad = function(y) -1/y ^ 2,
    y.hess = function(x) 2/(x ^ 3),
    x.hess = function(y) 2/(y ^ 3),
    x.in.domain = function(x) all(abs(x) > .Machine$double.eps),
    y.in.domain = function(y) all(abs(y) > .Machine$double.eps)
  )
}

#' Logit reparametrisation (1D)
#'
#' y = logit(x)
#'
#' @return logit function as \code{\link{ExpFam_dist_reparam}}
#' @export
Reparam_Logit <- function(){
  z <- ExpFam_dist_reparam(
    y.from.x = function(x) log(x) - log(1 - x),
    x.from.y = function(y) 1/(1 + exp(-y)),
    y.grad = function(x) 1/(x*(1 - x)),
    x.grad = function(y) (exp(y)/(1 + exp(y)))/(1 + exp(y)),
    y.hess = function(x) 1/((x - 1)^2) - 1/(x^2),
    x.hess = function(y) -exp(y)*(exp(y) - 1)/(1 + exp(y))^3,
    x.in.domain = function(x) x > 0 & x < 1,
    y.in.domain = function(y) TRUE
  )
  z$stats.glm.link <- stats::make.link("logit")
  z <- N_autogenerate_from_stat_GLM(z)

  z
}


