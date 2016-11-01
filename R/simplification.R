#' Function reparametrization
#'
#' Creates a function being a reparametrisation (with simplification) of
#'  \code{f1 * f2}
#'
#' @param f1 function 1
#' @param f2 function 2
#' @param res.proto resulting function prototype with NULL body, e.g. \code{function(theta) NULL}
#'
#' @examples
#'
#'  reparametrize_function_with_proto(
#'    function(x) x^2,
#'    function(y) y^3,
#'    function(z) NULL
#'    )
#'    # function(z) z^6
#'
#'  reparametrize_function_with_proto(
#'    function(x) x^2 + 2*x,
#'    function(y) (y+1)^3,
#'    function(z) NULL
#'    )
#'    # function(z) ((1 + z)^3 + 2) * (1 + z)^3
#'
#' @export
reparametrize_function_with_proto <- function(f1, f2, res.proto){

  f1s <- Deriv::Simplify(f1)
  f2s <- Deriv::Simplify(f2)
  f1s.param.sub <- list(quote(xxx.tmp.xxx.tmp))  # tmp name
  res.param.sub <- list(as.symbol(methods::formalArgs(res.proto)))  # res function param name
  #list(as.symbol("ljkfsda"))

  #tmp <- list(formals(res.proto)[1])

  names(f1s.param.sub) <- names(formals(f1s)[1])
  names(res.param.sub) <- names(formals(f2s)[1])
  f1s.q <- pryr::substitute_q(body(f1s), f1s.param.sub)
  f2s.q <- pryr::substitute_q(body(f2s), res.param.sub)
  f3s <- pryr::substitute_q(f1s.q, list(xxx.tmp.xxx.tmp = f2s.q))
  f3ss <- Deriv::Simplify(f3s)
  #pryr::make_function(list(x = NULL), body = quote(1/x^2))
  as.function(
    c(
      formals(res.proto),
      #as.list(as.symbol(res.param.name)),
      f3ss
    ),
    envir = environment(f1))
}

#' @export
rename_farg_with_proto <- function(f, res.proto){
  res.param.sub <- list(as.symbol(methods::formalArgs(res.proto)))
  names(res.param.sub) <- names(formals(f)[1])
  f.q <- pryr::substitute_q(body(f), res.param.sub)
  as.function(
    c(
      formals(res.proto),
      #as.list(as.symbol(res.param.name)),
      f.q
    ),
    envir = environment(f))

}

#' @export
reparametrize_function_grad_with_proto <- function(
  f.from.z, f.from.z.grad,
  z.from.x, z.from.x.grad,
  res.proto, simplify = TRUE
){
  if (simplify) {
    f.from.z.grad.sub <- reparametrize_function_with_proto(
      f1 = f.from.z.grad,
      f2 = z.from.x,
      res.proto = res.proto
    )

    z.from.x.grad.ren <- rename_farg_with_proto(
      f = z.from.x.grad,
      res.proto = res.proto
    )

    f5s.q <- pryr::subs(
      A*B,
      list(
        A = body(f.from.z.grad.sub),
        B = body(z.from.x.grad.ren)
      )
    )
    f5s <- Deriv::Simplify(f5s.q)
    ret <- as.function(
      c(
        formals(res.proto),
        f5s
      ),
      envir = environment(f.from.z.grad))

  } else {
    ret <- function(x) f.from.z.grad(z.from.x(x))*z.from.x.grad(x)
  }
  ret
}

#' @export
reparametrize_function_hess_with_proto <- function(
  f.from.z, f.from.z.grad, f.from.z.hess,
  z.from.x, z.from.x.grad, z.from.x.hess,
  res.proto, simplify = TRUE
){
  if (simplify) {
    z.from.x.grad.sq <- reparametrize_function_with_proto(
      f1 = function(x) x ^ 2,
      f2 = z.from.x.grad,
      res.proto = res.proto
    )

    f.hess <- reparametrize_function_with_proto(
      f1 = f.from.z.hess,
      f2 = z.from.x,
      res.proto = res.proto
    )

    f.from.z.grad.sub <- reparametrize_function_with_proto(
      f1 = f.from.z.grad,
      f2 = z.from.x,
      res.proto = res.proto
    )

    z.from.x.hess.ren <- rename_farg_with_proto(
      f = z.from.x.hess,
      res.proto = res.proto
    )


    f5s.q <- pryr::subs(
      A*B + C*D,
      list(
        A = body(z.from.x.grad.sq),
        B = body(f.hess),
        C = body(f.from.z.grad.sub),
        D = body(z.from.x.hess.ren)
      )
    )
    f5s <- Deriv::Simplify(f5s.q)
    ret <- as.function(
      c(
        formals(res.proto),
        f5s
      ),
      envir = environment(f.from.z.grad))

  } else {
    ret <- function(x){
      z.from.x.grad(x) ^ 2 * f.from.z.hess(z.from.x(x)) +
        f.from.z.grad(z.from.x(x)) * z.from.x.hess(x)
      }
  }
  ret
}
