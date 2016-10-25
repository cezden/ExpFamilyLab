#' @export
f_num_opt_proto <- function(num.opt, num.opt.TRUE, num.opt.FALSE, error.cond, err.msg){
  f.error <- if (error.cond) {
      function() stop(err.msg)
    } else {
      function_NULL
      }
  functor_get_first_not_null_preference_guarded_eval(
    num.opt,
    num.opt.TRUE,
    num.opt.FALSE,
    f.error
  )
}

#' @export
f_num_opt_generic <- function(obj, fun.name, num.opt = TRUE, error.cond = TRUE){
  fun.name.opt <- paste0("N.", fun.name)
  f_num_opt_proto(num.opt, c(obj[[fun.name.opt]]), c(obj[[fun.name]]), error.cond, paste("No", fun.name))
}

#' @export
f_eta_from_theta <- function(obj, num.opt = TRUE, error.cond = TRUE){
  f_num_opt_generic(obj, "eta.from.theta", num.opt, error.cond)
}

#' @export
f_theta_from_eta <- function(obj, num.opt = TRUE, error.cond = TRUE){
  f_num_opt_proto(num.opt, c(obj$N.theta.from.eta), c(obj$theta.from.eta), error.cond, "No theta.from.eta")
}

#' @export
f_B_from_theta <- function(obj, num.opt = TRUE, error.cond = TRUE){
  f_num_opt_proto(num.opt, c(obj$N.B.from.theta), c(obj$B.from.theta), error.cond, "No B.from.theta")
}

#' @export
f_A_from_eta <- function(obj, num.opt = TRUE, error.cond = TRUE){
  f_num_opt_proto(num.opt, c(obj$N.A.from.eta), c(obj$A.from.eta), error.cond, "No A.from.eta")
}

