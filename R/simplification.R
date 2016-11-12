#' Function reparametrization
#'
#' Creates a function being a reparametrisation (with simplification) of
#'  \code{f1 * g}, e.g. f1(g(x))
#'
#' @param f1 function 1
#' @param g function 2
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
reparametrize_function_with_proto <- function(f1, g, res.proto){

  #f1s <- Deriv::Simplify(f1)
  #f2s <- Deriv::Simplify(g)
  f1s <- multisimplify.expression(f1)
  f2s <- multisimplify.expression(g)
  f1s.param.sub <- list(quote(xxx.tmp.xxx.tmp))  # tmp name
  res.param.sub <- list(as.symbol(methods::formalArgs(res.proto)))  # res function param name
  #list(as.symbol("ljkfsda"))

  #tmp <- list(formals(res.proto)[1])

  names(f1s.param.sub) <- names(formals(f1s)[1])
  names(res.param.sub) <- names(formals(f2s)[1])
  f1s.q <- pryr::substitute_q(body(f1s), f1s.param.sub)
  f2s.q <- pryr::substitute_q(body(f2s), res.param.sub)
  f3s <- pryr::substitute_q(f1s.q, list(xxx.tmp.xxx.tmp = f2s.q))

  f3ss <- multisimplify.expression(f3s)
  #f3ss <- Deriv::Simplify(f3s)

  #pryr::make_function(list(x = NULL), body = quote(1/x^2))
  as.function(
    c(
      formals(res.proto),
      #as.list(as.symbol(res.param.name)),
      f3ss
    ),
    envir = environment(f1))
}


#' Change the name of 1st formal argument
#' @param f (function)
#' @param res.proto "prototype" function with desired formal name, e.g. \code{function(x) NULL}
#' @return function \code{f} with resubstituted first argument
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

#' Gradient of function reparametrization
#'
#' Function calculates the gradient df(z(x))/dx of f(z(x))
#'
#' @param f.from.z (\code{function}) f(z)
#' @param f.from.z.grad (\code{function}) df/dz (z)
#' @param z.from.x (\code{function}) z(x)
#' @param z.from.x.grad (\code{function}) dz/dx (x)
#' @param res.proto "prototype" function with desired formal name, e.g. \code{function(x) NULL}
#' @param simplify (\code{logical}) should the expression simplification be used?
#' @export
reparametrize_function_grad_with_proto <- function(
  f.from.z = NULL, f.from.z.grad,
  z.from.x, z.from.x.grad,
  res.proto, simplify = TRUE
){
  if (simplify) {
    f.from.z.grad.sub <- reparametrize_function_with_proto(
      f1 = f.from.z.grad,
      g = z.from.x,
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
    #f5s <- Deriv::Simplify(f5s.q)
    f5s <- multisimplify.expression(f5s.q)
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

#' Hessian of function reparametrization
#'
#' Function calculates the hessian d{^2}f(z(x))/dx^2 of f(z(x))
#'
#' @param f.from.z (\code{function}) f(z)
#' @param f.from.z.grad (\code{function}) df/dz (z)
#' @param f.from.z.hess (\code{function}) d^2 f/dz^2 (z)
#' @param z.from.x (\code{function}) z(x)
#' @param z.from.x.grad (\code{function}) dz/dx (x)
#' @param z.from.x.hess (\code{function}) d^2 z/dx^2 (x)
#' @param res.proto "prototype" function with desired formal name, e.g. \code{function(x) NULL}
#' @param simplify (\code{logical}) should the expression simplification be used?
#' @export
reparametrize_function_hess_with_proto <- function(
  f.from.z = NULL, f.from.z.grad, f.from.z.hess,
  z.from.x, z.from.x.grad, z.from.x.hess,
  res.proto, simplify = TRUE
){
  if (simplify) {
    z.from.x.grad.sq <- reparametrize_function_with_proto(
      f1 = function(x) x ^ 2,
      g = z.from.x.grad,
      res.proto = res.proto
    )

    f.hess <- reparametrize_function_with_proto(
      f1 = f.from.z.hess,
      g = z.from.x,
      res.proto = res.proto
    )

    f.from.z.grad.sub <- reparametrize_function_with_proto(
      f1 = f.from.z.grad,
      g = z.from.x,
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

    #f5s <- Deriv::Simplify(f5s.q)
    f5s <- multisimplify.expression(f5s.q)

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

#' @export
multisimplify.expression <- function(expr){
  expr.nchar <- nchar(paste(deparse(expr), collapse = " "))
  if (base::requireNamespace("Ryacas", quietly = TRUE)) {
    expr.tmp <- simplify_Ryacas(expr)
    expr.tmp.nchar <- nchar(paste(deparse(expr.tmp), collapse = " "))
    if (expr.nchar * 1.1 >= expr.tmp.nchar) {
      print("-----------Ryacas")
      expr <- expr.tmp
      expr.nchar <- expr.tmp.nchar
    } else {
      print(expr)
      print(expr.tmp)
      print("-----------")
    }
  } else {
    paste("No Ryacas...")
  }
  expr.tmp <- if (is.function(expr)) {
    Deriv::Simplify(expr, env = environment(expr))
  } else {
    Deriv::Simplify(expr)
  }
  expr.tmp.nchar <- nchar(paste(deparse(expr.tmp), collapse = " "))
  if (expr.nchar * 1.1 >= expr.tmp.nchar) {
    expr.tmp
  } else {
    expr
  }
}

#' Ryacas expression simplifyier
#'
#' Based on Deriv::Simplify
#' @export
simplify_Ryacas <- function(expr){

  if (!base::requireNamespace("Ryacas", quietly = TRUE)) {
    return(expr)
  }

  if (is.expression(expr)) {
    Ryacas::Simplify.default(expr) %>% Ryacas::as.expression.Sym()
  }
  else if (is.function(expr)) {
    as.function(
      c(
        formals(expr),
        as.list(Ryacas::Simplify.default(Ryacas::bodyAsExpression(expr)) %>% Ryacas::as.expression.Sym())
      ),
      envir = environment(expr)
    )
  } else if (is.call(expr)) {
    #warning("call object --> expression")
    expr.tmp <- Ryacas::Simplify.default(as.expression(expr)) %>% Ryacas::as.expression.Sym()
    expr.tmp[[1]]
  } else {
    expr
  }

}