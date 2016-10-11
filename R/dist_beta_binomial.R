distributions.beta <- function(alpha, beta){

}

#' The mode of beta distribution
distributions.mode.beta <- function(beta.coeffs){
  beta.mode.vals <- function(alpha, beta){
    if (is.na(alpha) | is.na(beta)){
      return (NA)
    }
    if ((alpha>1 && beta==1) || (alpha>=1 && beta<1)){
      return(1)
    }

    if ((alpha == 1 && beta > 1) || (alpha < 1 && beta >= 1)){
      return(0)
    }
    if (alpha<1 && beta<1){
      return(c(0,1))
    }
    mode.m <- ((alpha - 1) / (alpha + beta - 2))
    return(mode.m)
  }

  with(beta.coeffs, {
    ## alpha, beta could be vectors...
    return(mapply(beta.mode.vals, alpha, beta))
  })
}

distributions.mean.beta <- function(beta.coeffs){
  with(beta.coeffs,{
    return ((alpha) / (alpha+beta))
  })
}

function(){
  distributions.mode.beta(list(list(alpha=0.5, beta=0.3)))
  distributions.mode.beta(list(alpha=c(0.5, 0.4), beta=c(0.3, 2)))[[1]]
}