tmptmp <- function(){
  nn <- (function() {
    qq <- 1
    function(x) x*qq
  })()
  pryr::unenclose(nn)

  mm <- pryr::unenclose(z.mean$B.hess)

  body(mm)

  Deriv::Simplify(substitute(1/theta^2, list(theta = quote(1/x))))

  Deriv::Simplify("1/(1/(1/x)^2)")

  ls(Deriv::simplifications)


  ttt <- as.list(Deriv::simplifications)
  ttt$sqrt
  ttt$sqrt(as.call(list(sqrt, quote(sqrt(x)))))
  Deriv::Simplify("sqrt(sqrt(a+1))")
  simplr::simplifyq(sqrt(sqrt(a+1)))
  Deriv::Simplify("f(x^(-1)) + f(1/x)")
  simplr::simplifyq(f(x^(-1)) + f(1/x))

  Deriv::Simplify("log(2+x-1)")
  Deriv::Simplify("log(2^x)")
  exp

  list(
    c("(+ x x)", "(* 2 x)"),
    c("(+ 0 x)", "x"),

    c("(log (+ 1 x) )", "(log1p x)"),
    c("(log (^ x p) )", "(* (log p) x )"),
    c("(log (exp x) )", "(x)"),
    c("(log (exp x ) p )", "(/ x (log p))"),
    c("(log (sqrt x) )", "(* 0.5 (log x) )"),
    c("(log () )", "()"),
    c("(log () )", "()"),
    c("(log () )", "()"),
    c("(sqrt (sqrt x) )", "(^ x 0.25)"),
    c("(sqrt (* x x) )", "(abs x)"),
  )
  ??SEXP


  parse(text = deparse(quote(f(a,b,c,d)))) %>% getParseData
  pryr::ast(f(a,b,c,d))
  pryr::ast(f(a,b,c,d))
  codetools::showTree(quote(f(a,b,c,d)))
  codetools::showTree(quote(log(exp(x), 3)))

  codetools::makeCodeWalker()
  function (
    ...,
    handler = function(v, w) NULL,
    call = function(e, w)
      for (ee in as.list(e))
        if (!missing(ee))
          walkCode(ee, w),
    leaf = function(e, w) print(e)
  )
    list(handler = handler, call = call, leaf = leaf, ...)

  walkCode_ <- function (e, w = makeCodeWalker())
  {
    if (typeof(e) == "language") {
      if (typeof(e[[1]]) %in% c("symbol", "character")) {
        h <- w$handler(as.character(e[[1]]), w)
        if (!is.null(h))
          h(e, w)
        else w$call(e, w)
      }
      else w$call(e, w)
    }
    else w$leaf(e, w)
  }

  z <- expression(a+1)
  zz <- z[[1]]
  length(zz)
  zz[[1]] %>% typeof
  zz[[1]] %>% as.character
  zz[[2]] %>% typeof
  zz[[3]] %>% typeof
  str(zz)

  z <- expression(a+1+log(x^3,2)-min(x,y))

  zz <- z[[1]]
  length(zz)
  zz[[1]]
  zz[[2]]
  zz[[3]]

  zz[[1]] %>% typeof
  zz[[1]] %>% as.character
  zz[[2]] %>% typeof
  zz[[3]] %>% typeof
  str(zz)



  codetools::showTree()

  showTree_ <- function (e){
    w <- makeCodeWalker(call = showTreeCall, leaf = showTreeLeaf,
                        write = write)
    codetools::walkCode(e, w)
    w$write("\n")
  }

  showTreeLeaf_ <- function (e, w)
  {
    if (typeof(e) == "symbol") {
      if (e == "(")
        w$write("\"(\"")
      else if (e == "{")
        w$write("\"{\"")
      else w$write(e)
    }
    else w$write(deparse(e))
  }


  showTreeCall_ <-
    function (e, w)
    {
      w$write("(")
      walkCode(e[[1]], w)
      for (a in as.list(e[-1])) {
        w$write(" ")
        if (missing(a))
          w$write("<Missing>")
        else walkCode(a, w)
      }
      w$write(")")
    }






  ex3 <- expression(u, v, 1+ 0:9)
  ex3[[3]]




  class(pryr::ast(1/(1/x)^2))
  simplr::simplify(parse(text="1/(1/x)^2"))
  simplr::simplifyq(1/(1/(1/x)^2))
  simplr::simplifyq((1/x)^2)
  simplr::simplifyq(x + x + x)
  simplr::simplify(body(z.mean$B.hess))

  parse(text = deparse(quote(-x))) %>% getParseData


  codetools::walkCode(quote(1/(1/x)^2))
  codetools::showTree(quote(1/(1/x)^2))
  parse(text = deparse(quote(1/(1/x)^2))) %>% getParseData

  parse(text = deparse(quote(x+x+x))) %>% getParseData
  codetools::showTree(quote(x+x+x))
  codetools::showTree(quote(-x))
  pryr::ast(x+x+x)
  pryr::ast(-x)
  simplr::simplifyq(x + x + x)
  simplr::simplifyq(x^2 + 2*x + 1)
  Deriv::Simplify("x^2 + 2*x + 1")


  as.list(environment(nn))
  substitute(expression(nn(x)), environment(nn))
  Deriv::Simplify(nn, environment(nn))
  getAnywhere( Simplify_ )
  getAnywhere( Simplify.rule )

  ex1 <- enquote(body(z.mean$eta.from.theta))

  nnn <- body(pryr::unenclose(z.mean$B.hess))
  class(nnn)
  pryr::partial(z.mean$B.hess, theta = 1/x)

  class(ex1)
  substitute(-1/theta, list(theta = 1) )
  Deriv::Simplify(z.mean$eta.from.theta, environment(z.mean$eta.from.theta))
  Deriv::Simplify(z.mean.inv$B.hess, environment(z.mean.inv$B.hess))

  z.mean.inv$B.hess(2)

  formals(z.mean$eta.from.theta)
  environment(z.mean$eta.from.theta)

}

tmptmp_inverse <- function(){
  z <- ExpFam_dist_canonical(
    h = function(x) (x >= 0) * 1.0,
    A.from.eta = function(eta) -log(-eta), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1,
    A.grad = function(eta) -1/eta,
    A.hess = function(eta) 1/eta ^ 2,
    eta.in.domain = function(eta) all(eta < 0)
  )
  z.rate <- ExpFam_param(
    org.dist = z,
    eta.from.theta = function(theta) -theta,
    theta.from.eta = function(eta) -eta,
    B.from.theta = function(theta) -log(theta),
    B.grad = function(theta) -1/theta,
    B.hess = function(theta) 1/theta ^ 2,
    theta.in.domain = function(theta) all(theta > 0),
    param.type = "rate"
  )



  z <- z.mean
  z <- z.rate

  z.mean.inv <- ExpFam_param(
    org.dist = z.mean$org.dist,
    eta.from.theta = function(theta) z.mean$eta.from.theta(1/theta),
    theta.from.eta = function(eta) 1/z$theta.from.eta(eta),
    B.from.theta = function(theta) z$B.from.theta(1/theta),
    B.grad = function(theta) -z$B.grad(1/theta)/(theta^2),
    B.hess = function(theta) (2*theta*z.mean$B.grad(1/theta) + z.mean$B.hess(1/theta))/theta^4,
    theta.in.domain = function(theta) z$theta.in.domain(1/theta),
    param.type = "mean-inv"
  )

  z.t <- Reparam_Inverse()

  ExpFam_reparametrize.ExpFam_param(z.rate, z.t, "rate-inv")



  z <- z.rate
  z.mean.inv2 <- ExpFam_param(
    org.dist = z$org.dist,
    eta.from.theta = function(theta) z$eta.from.theta(z.t$y.from.x(theta)),
    theta.from.eta = function(eta) z.t$x.from.y(z$theta.from.eta(eta)),
    B.from.theta = function(theta) z$B.from.theta(z.t$y.from.x(theta)),
    B.grad = function(theta) z$B.grad(z.t$y.from.x(theta))*z.t$y.grad(theta),
    B.hess = function(theta) (2*theta*z$B.grad(1/theta) + z$B.hess(1/theta))/theta^4,
    theta.in.domain = function(theta) z.t$x.in.domain(theta) && z$theta.in.domain(z.t$y.from.x(theta)),
    param.type = "mean-inv"
  )

  f1s <- Deriv::Simplify(z$eta.from.theta)
  f1s.q <- pryr::substitute_q(body(f1s), list(theta = quote(xxx.tmp.xxx.tmp)))
  f2s <- Deriv::Simplify(z.t$y.from.x)
  f2s.q <- pryr::substitute_q(body(f2s), list(x = quote(theta)))
  f3s <- pryr::substitute_q(f1s.q, list(xxx.tmp.xxx.tmp = f2s.q))
  f3ss <- Deriv::Simplify(f3s)

  f1s <- Deriv::Simplify(z$B.grad)
  f1s.q <- pryr::substitute_q(body(f1s), list(theta = quote(xxx.tmp.xxx.tmp)))
  f2s <- Deriv::Simplify(z.t$y.from.x)
  f2s.q <- pryr::substitute_q(body(f2s), list(x = quote(theta)))
  f3s <- pryr::substitute_q(f1s.q, list(xxx.tmp.xxx.tmp = f2s.q))
  f3ss <- Deriv::Simplify(f3s)

  f4s <- Deriv::Simplify(z.t$y.grad)
  f4s.q <- pryr::substitute_q(body(f4s), list(x = quote(theta)))
  f5s.q <- pryr::subs(B.grad*y.grad, list(B.grad = f3ss, y.grad = f4s.q))
  f5ss <- Deriv::Simplify(f5s.q)


  str(quote(theta))
  qqz <- "theta"
  as.symbol(qqz)

  reparametrize_function_1 <- function(f1, f2, res.param.name){
    f1s <- Deriv::Simplify(f1)
    f2s <- Deriv::Simplify(f2)
    f1s.param.sub <- list(quote(xxx.tmp.xxx.tmp))  # tmp name
    res.param.sub <- list(as.symbol(res.param.name))  # res function param name

    names(f1s.param.sub) <- names(formals(f1s)[1])
    names(res.param.sub) <- names(formals(f2s)[1])
    f1s.q <- pryr::substitute_q(body(f1s), f1s.param.sub)
    f2s.q <- pryr::substitute_q(body(f2s), res.param.sub)
    f3s <- pryr::substitute_q(f1s.q, list(xxx.tmp.xxx.tmp = f2s.q))
    f3ss <- Deriv::Simplify(f3s)
    #pryr::make_function(list(x = NULL), body = quote(1/x^2))
    par.names <- list(xxx = NULL)
    names(par.names) <- res.param.name
    as.function(
      c(
        par.names,
        #as.list(as.symbol(res.param.name)),
        f3ss
      ),
      envir = environment(f1))
  }


  res.proto <- function(theta) NULL

  tststs <- reparametrize_function_1(function(x) x^2, z.t$y.grad, "theta")
  tststs_ <- reparametrize_function_2(function(x) x^2, z.t$y.grad, function(theta) NULL)
  tststs2 <- reparametrize_function_1(z$B.hess, z.t$y.from.x, "theta")

  tststs3 <- reparametrize_function_1(z$B.grad, z.t$y.from.x, "theta")
  tststs4 <- reparametrize_function_1(z.t$y.hess, function(x) x, "theta")
  f5s.q <- pryr::subs(
    A*B + C*D,
    list(
      A = body(tststs),
      B = body(tststs2),
      C = body(tststs3),
      D = body(tststs4)
    )
  )
  f5s <- Deriv::Simplify(f5s.q)

  simplr::simplify(tststs)
  getAnywhere(simplr::simplify)
  simplr::simplify(tststs2)
  simplr::simplify(tststs3)
  simplr::simplify(tststs4)

  Deriv::Simplify("(-(1/theta^2))^2/(1/theta)^2 - 2/theta^2")

  Deriv::Simplify("GBIS* (-1/x) + GPRSQ* (1/x^2)")
  Deriv::Simplify("(GPRSQ/x - GBIS)/x")

  simplr::simplifyq((-(theta^(-2)))^2/(theta)^(-2) - 2*theta^(-2))
  simplr::simplifyq((-theta^(-2))^2/(theta)^(-2) - 2*theta^(-2))
  Deriv::Simplify("(-theta^(-2))^2/(theta)^(-2) - 2*theta^(-2)")

  deparse(quote((-(theta^(-2)))^2/(theta)^(-2) - 2*theta^(-2)))

  ( -theta^(-2) )^2 / (theta)^(-2) - 2*theta^(-2)

  simplr::simplifyq(((theta^(-2)))^2/(theta)^(-2) - 2*theta^(-2))



  f4s <- Deriv::Simplify(z.t$y.grad)
  f4s.q <- pryr::substitute_q(body(f4s), list(x = quote(theta)))
  f5s.q <- pryr::subs(B.grad*y.grad, list(B.grad = f3ss, y.grad = f4s.q))
  f5ss <- Deriv::Simplify(f5s.q)




  z.t$y.grad(theta)




  f1s <- Deriv::Simplify(z$B.from.theta)
  f1s.q <- pryr::substitute_q(body(f1s), list(theta = quote(xxx.tmp.xxx.tmp)))
  f2s <- Deriv::Simplify(z.t$y.from.x)
  f2s.q <- pryr::substitute_q(body(f2s), list(x = quote(theta)))
  f3s <- pryr::substitute_q(f1s.q, list(xxx.tmp.xxx.tmp = f2s.q))



  f2s.b <- body(f2s) #getting 'language'
  is.expression(f2s.b)
  is.call(f2s.b)
  lapply(as.list(f2s.b)[-1], typeof)


  substitute(, list(x = quote(theta)))

  ff0 <- function(theta) z$eta.from.theta(theta = z.t$y.from.x(x = theta))

  f0 <- z.mean.inv2$eta.from.theta
  f0s <- Deriv::Simplify(ff0)
  f1 <- z$eta.from.theta
  f2 <- z.t$y.from.x
  f1s <- Deriv::Simplify(f1)
  f2s <- Deriv::Simplify(f2)
  Deriv::Simplify(quote(function(x) f1s(f2s(x))))
  Lincomb
}


tmptmp.yacas <- function(){

  library(Ryacas)
  ??Ryacas
  x <- Ryacas::Sym("x")
  sym0 <- x^3
  sym1 <- x^2+2*x^2
  sym2 <- 2 * sym0
  sym3 <- Sym(6) * pi * x
  sym4 <- sym1 * (1 - sin(sym3)) / sym2
  print(sym4)
  sym5 <- Ryacas::Simplify(sym4)
  print(sym5)
  class(sym5)
  class(as.expression(log(sym5)))
  qqq <- as.expression(log(sym5))
  eval(qqq, list(x=1:5/10))
  qqq <- Ryacas::bodyAsExpression(z.mean$B.from.theta)
  qqq <- Ryacas::bodyAsExpression(z.mean$B.grad)
  qqq <- Ryacas::bodyAsExpression(function(theta) log(1) - log(1-theta))

  qqq
  Ryacas::deriv.Expr(qqq, name = Expr(as.name("theta")))
  Ryacas::deriv.Expr(qqq, theta)
  Ryacas::deriv.Expr(log(1) - log(1-theta), theta)
  Ryacas::Simplify(qqq)


  qqq <- expression((-(1/theta^2))^2/(1/theta)^2 - 2/theta^2)
  expr.y <- Ryacas::Expr(qqq)
  expr.tmp <- Ryacas::Simplify.Expr(expr.y) %>% Ryacas::as.expression.Sym()

  multisimplify.expression(qqq)
  simplify_Ryacas(qqq)

  simplify_Ryacas(function(theta) (-(1/theta^2))^2/(1/theta)^2 - 2/theta^2)
  multisimplify.expression(function(theta) (-(1/theta^2))^2/(1/theta)^2 - 2/theta^2)

  fff <- function(theta) (-(1/theta^2))^2/(1/theta)^2 - 2/theta^2

  simplify_Ryacas(fff)
  multisimplify.expression(fff)

  simplify_Ryacas(qqq)
  multisimplify.expression(qqq)

  simplify_Ryacas(expression( log(exp(log(theta) - log(1 - theta)) + 1)))

  qqq <- expression( log(exp(log(theta) - log(1 - theta)) + 1))
  expr.y <- Ryacas::Expr(qqq)
  expr.tmp <- Ryacas::Simplify.Expr(expr.y) %>% Ryacas::as.expression.Sym()
  qqq <- function (theta)
    (-(1/theta^2))^2/(1/theta)^2 - 2/theta^2
  simplify_Ryacas(qqq)
  multisimplify.expression(qqq)
  as.call(qqq)


  ppp <- c(
    formals(fff),
    as.list(Ryacas::Simplify.default(Ryacas::bodyAsExpression(fff)) %>% Ryacas::as.expression.Sym())
  )
  str(ppp)
  class(ccc)
  Expr(qqq)

  as.function(ppp, envir = base::environment(fff))

  Ryacas::Simplify(qqq)
  qqq <- Ryacas::bodyAsExpression(function(theta) (-(1/theta^2))^2/(1/theta)^2 - 2/theta^2)
  class(qqq)
  qqqE <- Ryacas::Expr(qqq)
  qqq1 <- Ryacas::Simplify(qqq) %>% Ryacas::as.expression.Sym()
  qqq1E <- Ryacas::Simplify(qqqE) %>% Ryacas::as.expression.Sym()
  str(qqq1)
  qqq2 <- Deriv::Simplify(qqq1)
  qqq3 <- Ryacas::Simplify(qqq2) %>% Ryacas::as.expression.Sym()

  qqq <- as.expression(log(sym5))

  sss <- Ryacas::deriv.Sym()

  sss <- expression((-(1/theta^2))^2/(1/theta)^2 - 2/theta^2)
  expr.y <- Ryacas::Expr(sss)
  expr.tmp <- Ryacas::deriv.Sym(expr.y, name = "theta") %>% Ryacas::as.expression.Sym()
  str(expr.tmp)

}