tmptmp <- function(){
  getAnywhere("logLik.glm")
  obj.binomial <- binomial()
  str(obj.binomial)
  binomial(link = "logit")
  obj.gaussian <- gaussian(link = "identity")
  str(obj.gaussian)
  obj.gaussian$variance
  obj.gaussian$aic
  obj.gaussian$dev.resids
  obj.gaussian$mu.eta
  Gamma(link = "inverse")
  inverse.gaussian(link = "1/mu^2")
  poisson(link = "log")
  quasi(link = "identity", variance = "constant")
  quasibinomial(link = "logit")
  quasipoisson(link = "log")

  dbinom(1, 1, 0.5, log = TRUE) + dbinom(1-1, 1, 0.5, log = TRUE)
  dbinom(1, 1, 0.4, log = TRUE) + dbinom(1-1, 1, 0.6, log = TRUE)
  ttt.ttt <- function(y, p) y * log(p) + (1-y) * log(1 - p)
  ttt.ttt(1, 0.5)
  ttt.ttt(1, 0.00000001)/dbinom(1, 1, 0.00000001, log = TRUE)
  ttt.ttt(1, 1-0.0000001/2)/dbinom(1, 1, 1-0.0000001/2, log = TRUE)

  context("Distributions - Bernoulli")

  test_that("basic tests",{
    expect_is(distributions.Bernoulli(), "distributions.exp.dist")

    expect_true(distributions.density.theta(distributions.Bernoulli(), 0.6)(1)==0.6)   #ok
    expect_true(distributions.density.eta(distributions.Bernoulli(), 0)(1)==0.5)   #ok

    theta.test <- 0.6
    d.t.06 <- distributions.density.theta(distributions.Bernoulli(), theta.test)
    expect_true(all(c(d.t.06(0) == 1-theta.test, d.t.06(1) == theta.test)))

    d.t.07 <- distributions.density.theta(distributions.Bernoulli(), 0.7)
    expect_true(all(c(d.t.07(0) == 1-0.7, d.t.07(1) == 0.7)))

    eta.test <- 0.1
    eta.test.theta <- distributions.Bernoulli()$theta.from.eta(eta.test)
    d.e.0 <- distributions.density.eta(distributions.Bernoulli(), eta.test)
    expect_true(sum((c(d.e.0(0) - (1-eta.test.theta), d.e.0(1) - eta.test.theta))^2) < 1e-16)


  })



}