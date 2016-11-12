# canonical form, natural parameter, and natural parameter space:
## (general form of one-parameter) exponential family: exp(eta(theta)*S(x) - B(theta))*h(x) [DasGupta11a, def 18.1]
### S(x) is called natural suff. stat. for the family [DasGupta11a, def 18.2]
## canonical one-parameter exponential family: exp(eta*S(x) - A(eta))*h(x) [DasGupta11a, def 18.3]
### eta is called natural parameter

## natural exponential family (NEF): f(x|eta) = exp(eta'x - A(eta))*h(x)

###
# log1p(x) computes log(1+x) accurately also for |x| << 1.
# expm1(x) computes exp(x) - 1 accurately also for |x| << 1.


#' @export
ExpFam_dist_ext_Bernoulli <- function(){
  ## natural parameter (eta) logit(x)
  z <- ExpFam_dist_canonical(
    h = function(x) 1*(x == 0 | x == 1), #measure
    A.from.eta = function(eta) log(1 + exp(eta)), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1,
    A.grad = function(eta) 1/(1 + exp(-eta)),
    A.hess = function(eta) (exp(eta)/(1 + exp(eta)))/(1 + exp(eta)),
    eta.in.domain = function(eta) TRUE,
    family.name = "Bernoulli"
  )
  ## parameter "by convention": probability of success (p)
  ## theta: mean parametrisation == p

  z.mean <- ExpFam_reparametrize(
    z,
    Reparam_Logit(),
    hints = list(
      B.from.theta = function(theta) -log(1 - theta), #log-normaliser for mean parametrisation
      B.grad = function(theta) 1 / (1 - theta),
      B.hess = function(theta) 1 / (-1 + theta) ^ 2
      #
      #B.grad = function(theta) -1 / (1 + theta),
      #B.hess = function(theta) 1 / (1 + theta) ^ 2
    ),
    param.type = "mean"
  )

  # z.mean <- ExpFam_param(
  #   org.dist = z,
  #   eta.from.theta = function(theta) log(theta) - log(1 - theta),
  #   theta.from.eta = function(eta) 1/(1 + exp(-eta)),
  #   B.from.theta = function(theta) -log(1 - theta), #log-normaliser for mean parametrisation
  #   B.grad = function(theta) -1 / (1 + theta),
  #   B.hess = function(theta) 1 / (1 + theta) ^ 2,
  #   theta.in.domain = function(theta) all(theta > 0 & theta < 1),
  #   param.type = "mean"
  # )

  z.mean$stats.glm <- binomial()

  z$conjugate.prior.elems_raw <- list(
    h = function(x) 1,
    #A.from.eta = function(eta) lgamma(eta[, 2] - eta[, 1]) + lgamma(eta[, 1]) - lgamma(eta[, 2]),
    A.from.eta = function(eta) lbeta(eta[, 2] - eta[, 1], eta[, 1]),
    A.grad = function(eta) {
      zz <- digamma(eta[, 2] - eta[, 1])
      cbind(digamma(eta[, 1]) - zz, zz - digamma(eta[, 2]))
    },
    A.hess = function(eta) {
      zz <- trigamma(eta[2] - eta[1])
      matrix(
        c(
          zz + trigamma(eta[1]), -zz,
          -zz, zz - trigamma(eta[2])
        ),
        nrow = 2,
        byrow = TRUE
      )
    },
    eta.in.domain = function(eta) all(eta[, 2] > eta[, 1] & eta[, 1] > 0)
  )

  ExpFam_dist_ext(
    canonical.param = z,
    parametrisations = list(z.mean)
  )
}

#' Bernoulli distribution
#'
#' @return object of class \code{\link{ExpFam_dist_ext}}
#' @export
dist_ext_Bernoulli <- function(){
  z <- ExpFam_dist_ext_Bernoulli()
  z[["mean"]] <- N_autogenerate_from_stat_GLM(z[["mean"]])

  ftheta <- f_theta_from_eta(z[["mean"]])
  z[["canonical"]]$N.density.eta <-
    function(x, eta, log = FALSE)
      dbinom(x = x, size = 1, prob = ftheta(eta), log = log)

  z[["mean"]]$N.density.theta <-
    function(x, theta, log = FALSE)
      dbinom(x = x, size = 1, prob = theta, log = log)

  #z <- N_autogenerate_functions(z)
  #class(z) <- c(class(z), "N_ExpFam_dist")
  z
}

#' @export
ExpFam_dist_ext_Exponential <- function(){
  # 'minus rate'
  z <- ExpFam_dist_canonical(
    h = function(x) (x >= 0) * 1.0,
    A.from.eta = function(eta) -log(-eta), #log-normaliser for natural parametrisation
    S = function(x) x, # Suff. stats
    eta.dim = 1,
    A.grad = function(eta) -1/eta,
    A.hess = function(eta) 1/eta ^ 2,
    eta.in.domain = function(eta) all(eta < 0),
    family.name = "Exponential"
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
  # z.mean <- ExpFam_param(
  #   org.dist = z,
  #   eta.from.theta = function(theta) -1/theta,
  #   theta.from.eta = function(eta) -1/eta,
  #   B.from.theta = function(theta) log(theta),
  #   B.grad = function(theta) 1/theta,
  #   B.hess = function(theta) -1/theta ^ 2,
  #   theta.in.domain = function(theta) all(theta > 0),
  #   param.type = "mean"
  # )
  z.rate.inv <- ExpFam_reparametrize(z.rate, Reparam_Inverse(), "mean")

  ExpFam_dist_ext(
    canonical.param = z,
    parametrisations = list(z.rate, z.rate.inv)
  )
}

#' Exponential distribution
#'
#' @return object of class \code{\link{ExpFam_dist_ext}}
#' @export
dist_ext_Exponential <- function(){
  z <- ExpFam_dist_ext_Exponential()

  z[["rate"]]$N.density.theta <-
    function(x, theta, log = FALSE)
      dexp(x = x, rate = theta, log = log)

  z[["mean"]]$N.density.theta <-
    function(x, theta, log = FALSE)
      dexp(x = x, rate = 1/theta, log = log)


  z[["canonical"]]$N.density.eta <-
    function(x, eta, log = FALSE)
      dexp(x = x, rate = -eta, log = log)

  z
}

