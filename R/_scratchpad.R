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
  mgcv::negbin()
  mgcv::nb()

  mgcv::betar()
  mgcv::family.mgcv
  mgcv::fix.family.link(fam)
  fix.family.var(fam)
  fix.family.ls(fam)
  fix.family.qf(fam)
  fix.family.rd(fam)

  rstanarm::priors
  rstanarm::stan_glm


  dbinom(1, 1, 0.5, log = TRUE) + dbinom(1-1, 1, 0.5, log = TRUE)
  dbinom(1, 1, 0.4, log = TRUE) + dbinom(1-1, 1, 0.6, log = TRUE)
  ttt.ttt <- function(y, p) y * log(p) + (1-y) * log(1 - p)
  ttt.ttt(1, 0.5)
  ttt.ttt(1, 0.00000001)/dbinom(1, 1, 0.00000001, log = TRUE)
  ttt.ttt(1, 1-0.0000001/2)/dbinom(1, 1, 1-0.0000001/2, log = TRUE)


}