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
    S = function(x) x # Suff. stats
  )
  z$stats.glm <- binomial()
  z
}

#' @export
dist_Bernoulli <- function(){
  z <- ExpFam_dist_Bernoulli()
  # ????
  z$logLik <- function(x, prob) dbinom(x, 1, prob, log = TRUE)
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
    S = function(x) x # Suff. stats
  )
  z$stats.glm <- poisson()
  z
}
