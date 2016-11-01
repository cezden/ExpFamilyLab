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
