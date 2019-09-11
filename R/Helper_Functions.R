# Matthew Smith, 9/10/19
# In this file, we create all of the helper functions that we will need in this package


#' Gamma Helper Function
#' 
#' @description This function represents the gamma function in my derivation, NOT the gamma function related to the factorial function.
#' @param x A real number (can be vectorized)
#' @return A real number (or vector) equalling (1-x/sqrt(2)) / (1+x/sqrt(2))
#' @export
g_helper <- function(x) {(1 - x/sqrt(2)) / (1 + x/sqrt(2))}


#' Beta Helper Function
#' 
#' @description This function represents the beta function in my derivation, NOT the beta function defined as an integral.
#' @param x A real number
#' @param y A real number
#' @return A real number equal to sqrt(2)*x/y
#' @export
b_helper <- function(x, y) {sqrt(2)*x/y}


#' Joint CDF Under Asymmetric Exponential Assumption
#' 
#' @description This function is the joint CDF of Delta and Zeta assuming the underlying data X_1 and X_2 follow independent exponential distributions with different parameters.
#' @param d The aboslute difference between two replicates
#' @param z The coefficient of variation between two replicates
#' @param lam_1 The (positive) parameter for the first replicate's exponential distribution
#' @param lam_2 The (positive) parameter for the second replicate's exponential distribution
#' @return The joint probability of observing Delta <= d and Z <= z
#' @export
F_exp_asym <- function(d, z, lam_1, lam_2) { #Joint CDF with asymmetric exponential assumption
  if(z < 0){
    return(0)
  }else if(d <= 0 && z >= 0 && z <= sqrt(2)){
    return(
      lam_1*lam_2*(1-g_helper(z))*exp((lam_1*g_helper(z)+lam_2)*d/(1-g_helper(z))) /
        ((lam_1*g_helper(z)+lam_2)*(lam_1+lam_2))
    )
  }else if(d <= 0 && z > sqrt(2)){
    return(lam_1*exp(lam_2*d)/(lam_1+lam_2))
  }else if(d > 0 && z >= 0 && z <= sqrt(2)){
    return(
      -2*lam_1*lam_2*exp(-((1+sqrt(2)/z)*lam_1+(sqrt(2)/z-1)*lam_2)*d/2)/(lam_1^2-lam_2^2) -
        lam_1*lam_2*(1+g_helper(z))*exp(-(lam_1+lam_2*g_helper(z))*sqrt(2)*d/(z*(g_helper(z)+1)))/
        ((lam_1+lam_2*g_helper(z))*(lam_2-lam_1)) +
        lam_1*lam_2*(1-(g_helper(z))^2)/((lam_1+lam_2*g_helper(z))*(lam_1*g_helper(z)+lam_2))
    )
  }else if(d > 0 && z > sqrt(2)){
    return(1 - lam_2*exp(-lam_1*d)/(lam_1+lam_2))
  }
}























































