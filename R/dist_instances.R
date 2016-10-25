# canonical form, natural parameter, and natural parameter space:
## (general form of one-parameter) exponential family: exp(eta(theta)*S(x) - B(theta))*h(x) [DasGupta11a, def 18.1]
### S(x) is called natural suff. stat. for the family [DasGupta11a, def 18.2]
## canonical one-parameter exponential family: exp(eta*S(x) - A(eta))*h(x) [DasGupta11a, def 18.3]
### eta is called natural parameter

## natural exponential family (NEF): f(x|eta) = exp(eta'x - A(eta))*h(x)

###



#' @export
ExpFam_dist_Bernoulli <- function(){
  ## parameter "by convention": probability of success (p)
  ## theta: mean parametrisation == p
  ## natural parameter (eta) logit(x)
  z <- ExpFam_dist(
    h = function(x) 1*(x == 0 | x == 1), #measure
    eta.from.theta = function(theta) log(theta) - log(1 - theta),
    theta.from.eta = function(eta) 1/(1 + exp(-eta)),
    B.from.theta = function(theta) -log(1 - theta), #log-normaliser for mean parametrisation
    A.from.eta = function(eta) log(1 + exp(eta)), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1
  )
  z$conjugate.prior.elems_eta_raw <- function() {
    list(
      h.c = function(x) 1,
      #A.from.eta = function(eta) lgamma(eta[, 2] - eta[, 1]) + lgamma(eta[, 1]) - lgamma(eta[, 2]),
      A.from.eta = function(eta) lbeta(eta[, 2] - eta[, 1], eta[, 1]),
      eta.valid = function(eta) all(eta[, 2] > eta[, 1] & eta[, 1] > 0)
    )
  }
  z$deta.from.theta_dtheta <- function(theta) 1/(theta*(1 - theta))
  z$dA.deta <- function(eta) exp(eta) / (1 + exp(eta)) # exp. val. of S
  z$d2A.deta2 <- function(eta) exp(eta) / (1 + exp(eta)) ^ 2 #exp. val. of S^2
  z$theta.from.params <- function(prob) prob
  z$stats.glm <- binomial()
  z
}

#' @export
dist_Bernoulli <- function(){
  z <- ExpFam_dist_Bernoulli()
  z <- N_autogenerate_from_stat_GLM(z)
  z$N.density.eta <- function(x, eta, log = FALSE) dbinom(x = x, size = 1, prob = z$N.theta.from.eta(eta), log = log)
  z$N.density.theta <- function(x, theta, log = FALSE) dbinom(x = x, size = 1, prob = theta, log = log)
  z <- N_autogenerate_functions(z)
  class(z) <- c(class(z), "N_ExpFam_dist")
  z
}


#' @export
ExpFam_dist_Poisson <- function(){
  ## parameter "by convention": lambda
  ## mean: lambda
  ## theta - lambda
  ## natural parameter (eta)
  z <- ExpFam_dist(
    h = function(x) (x == round(x))*(x >= 0)/gamma(x + 1), #measure  1/(x!)
    eta.from.theta = function(theta) log(theta),
    theta.from.eta = function(eta) exp(eta),
    B.from.theta = function(theta) theta, #log-normaliser for mean parametrisation
    A.from.eta = function(eta) exp(eta), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1
  )
  z$deta.from.theta_dtheta <- function(theta) 1/theta
  z$dA.deta <- function(eta) exp(eta)
  z$d2A.deta2 <- function(eta) exp(eta)
  z$theta.from.params <- function(lambda) lambda
  z$stats.glm <- poisson()
  z
}

#' @export
dist_Poisson <- function(){
  z <- ExpFam_dist_Poisson()
  z <- N_autogenerate_from_stat_GLM(z)
  #z$N.density.eta <- function(x, eta, log = FALSE) dpois(x, lambda = z$N.theta.from.eta(eta), log = log)
  z$N.density.theta <- function(x, theta, log = FALSE) dpois(x, lambda = theta, log = log)
  z <- N_autogenerate_functions(z)
  class(z) <- c(class(z), "N_ExpFam_dist")
  z
}


ExpFam_dist_Exponential_mean <- function(){
  ## parameterised by mean (mu)
  ## theta == mu:
  z <- ExpFam_dist(
    h = function(x) (x >= 0) * 1.0,
    eta.from.theta = function(theta) -1/theta,
    theta.from.eta = function(eta) -1/eta,
    B.from.theta = function(theta) -log(1/theta),
    A.from.eta = function(eta) -log(-eta), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1
  )

  z$eta.from.theta <- function(theta) -1/theta
  z$theta.from.eta <- function(eta) -1/eta

  z$grad.eta.from.theta <- function(theta) 1/(theta) ^ 2
  z$dA.deta <- function(eta) -1/eta # exp. val. of S
  z$d2A.deta2 <- function(eta) 1/eta ^ 2 #exp. val. of S^2
  z
}

#' @export
ExpFam_dist_Exponential_rate <- function(){
  ## parameter "by convention": rate (lambda)
  ## mean: 1/lambda
  ## theta - lambda
  ## natural parameter (eta)
  z <- ExpFam_dist(
    h = function(x) (x >= 0) * 1.0,
    eta.from.theta = function(theta) -theta,
    theta.from.eta = function(eta) -eta,
    B.from.theta = function(theta) -log(theta),
    A.from.eta = function(eta) -log(-eta), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1
  )
  z$dA.deta <- function(eta) -1/eta
  z$d2A.deta2 <- function(eta) 1/(eta ^ 2)
  z$theta.from.params <- function(lambda) lambda
  z$theta.from.mean <- function(mu) 1/mu
  #??vglmff
  #??gam
  z
}

#' @export
dist_Exponential <- function(){
  z <- ExpFam_dist_Exponential()
  #z$N.density.eta <- function(x, eta, log = FALSE) dexp(x, rate = z$theta.from.eta(eta), log = log)
  z$N.density.theta <- function(x, theta, log = FALSE) dexp(x, rate = theta, log = log)
  z <- N_autogenerate_functions(z)
  class(z) <- c(class(z), "N_ExpFam_dist")
  z
}

#' @export
ExpFam_dist_Beta <- function(){
  ## natural parameter (eta)

  #tmp <- mgcv::betar()
  #str(tmp)
  #tmp$linkfun(0.9)


  #non-canonical representation since S != x
  z <- ExpFam_dist(
    h = function(x) ifelse(x >= 0 & x <= 1, 1.0 / (x * (1 - x)), 0),
    eta.from.theta = NULL,
    theta.from.eta = NULL,
    B.from.theta = NULL,
    A.from.eta = function(eta) lbeta(eta[, 1], eta[, 2]),
    S = function(x) cbind(log(x), log(1 - x)), # Suff. stats
    eta.dim = 2
  )
  #z$dA.deta <- function(eta) -1/eta
  #z$d2A.deta2 <- function(eta) 1/(eta ^ 2)
  z$theta.from.params <- function(lambda) lambda
  z$theta.from.mean <- function(mu) 1/mu
  #??vglmff
  #??gam
  z
}


#' @export
dist_Beta <- function(){
  z <- ExpFam_dist_Beta()
  #z <- N_autogenerate_from_stat_GLM(z)
  z$N.density.eta <- function(x, eta, log = FALSE) dbeta(x = x, shape1 = eta[1], shape2 = eta[2], log = log)
  #z$N.density.theta <- function(x, theta, log = FALSE) dbinom(x = x, size = 1, prob = theta, log = log)
  #z <- N_autogenerate_functions(z)
  class(z) <- c(class(z), "N_ExpFam_dist")
  z
}


