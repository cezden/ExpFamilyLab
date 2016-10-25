#' @export
functor_eval_1 <- function(f, default.val, ...) {
  if (is.null(f)) {
    return(default.val)
  }
  f(...)
}

#' @export
functor_eval_first_not_null <- function(f.vec, default.val, ...) {
  f <- functor_get_first_not_null(f.vec, default.func = NULL)
  functor_eval_1(f, default.val = default.val, ...)
}


#' @export
function_NULL <- function(...) {
  NULL
}

#' @export
functor_NULL <- function(...) {
  function_NULL
}


#' @export
functor_get_first_not_null <- function(f.vec, default.func = function_NULL) {
  for (f in f.vec) {
    if (!is.null(f)) {
      return(f)
    }
  }
  return(default.func)
}

#' @export
functor_get_first_not_null_guarded_eval <- function(f.vec, default.func.eval = functor_NULL) {
  for (f in f.vec) {
    if (!is.null(f)) {
      return(f)
    }
  }
  return(default.func.eval())
}

#' @export
functor_get_first_not_null_preference_guarded_eval <- function(prefer.true, f.vec.true, f.vec.false, default.func.eval = functor_NULL) {
  f.vec <- if (prefer.true) {
    c(f.vec.true, f.vec.false)
  } else {
    f.vec.false
  }
  functor_get_first_not_null_guarded_eval(f.vec, default.func.eval = default.func.eval)
}

