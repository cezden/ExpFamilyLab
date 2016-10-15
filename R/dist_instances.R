#' @export
ExpFam_dist_Bernoulli <- function(){
  ## parameter "by convention": probability of success (p)
  ## theta: mean parametrisation == p
  ## natural parameter (eta) logit
  z <- ExpFam_dist(
    h = function(x) 1*(x == 0 | x == 1), #measure
    eta.from.theta = function(theta) log(theta) - log(1 - theta),
    theta.from.eta = function(eta) 1/(1 + exp(-eta)),
    B.from.theta = function(theta) -log(1 - theta), #log-normaliser for mean parametrisation
    A.from.eta = function(eta) log(1 + exp(eta)), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1
  )
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

#' @export
ExpFam_dist_Exponential <- function(){
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


