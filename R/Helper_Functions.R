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
F_exp_asym <- function(d, z, lam_1, lam_2) {
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


#' Joint CDF Under Symmetric Exponential Assumption
#' 
#' @description This function is the joint CDF of Delta and Zeta assuming the underlying data X_1 and X_2 follow independent and identically distributed exponential distributions.
#' @param d The aboslute difference between two replicates
#' @param z The coefficient of variation between two replicates
#' @param lam The (positive) parameter for the replicates's exponential distributions
#' @return The joint probability of observing Delta <= d and Z <= z
#' @export
F_exp_sym <- function(d, z, lam) {
  if(z < 0){
    return(0)
  }else if(d <= 0 && z >= 0 && z <= sqrt(2)){
    return((1-g_helper(z))*exp(lam*(1+g_helper(z))*d/(1-g_helper(z)))/(2*(1+g_helper(z))))
  }else if(d <= 0 && z > sqrt(2)){
    return(exp(lam*d)/2)
  }else if(d > 0 && z >= 0 && z <= sqrt(2)){
    return(
      lam*exp(-lam*b_helper(d,z))*d/2 +
        (1/2-1/(g_helper(z)+1))*(-1+exp(-lam*b_helper(d,z))+lam*b_helper(d,z)*exp(-lam*b_helper(d,z))) +
        (1-g_helper(z))/(2*(1+g_helper(z)))
    )
  }else if(d > 0 && z > sqrt(2)){
    return(1-exp(-lam*d)/2)
  }
}


#' Marginal CDF for Delta Under Asymmetric Exponential Assumption
#' 
#' @description This function is the marginal CDF of Delta (which follows an Asymmetric Laplace Distribution) assuming the underlying data X_1 and X_2 follow independent exponential distributions with different parameters.
#' @param d The absolute difference between two replicates
#' @param lam_1 The (positive) parameter for the first replicate's exponential distribution
#' @param lam_2 The (positive) parameter for the second replicate's exponential distribution
#' @return The marginal probability of observing Delta <= d
#' @export
F_D_asym <- function(d, lam_1, lam_2) {
  if(d < 0){
    return(
      (lam_1/(lam_1+lam_2)) * exp(lam_2*d)
    )
  }else{
    return(
      1 - (lam_2/(lam_1+lam_2))*exp(-lam_1*d)
    )
  }
}


#' Marginal CDF for Delta Under Symmetric Exponential Assumption
#' 
#' @description This function is the marginal CDF of Delta (which follows a Symmetric Laplace Distribution) assuming the underlying data X_1 and X_2 follow independent and identically distributed exponential distributions.
#' @param d The absolute difference between two replicates
#' @param lam The (positive) parameter for the replicates's exponential distributions
#' @return The marginal probability of observing Delta <= d
#' @export
F_D_sym <- function(d, lam) {
  if(d < 0){
    return(
      (1/2) * exp(lam*d)
    )
  }else{
    return(
      1 - (1/2)*exp(-lam*d)
    )
  }
}


#' Marginal CDF for Zeta Under Asymmetric Exponential Assumption
#' 
#' @description This function is the marginal CDF of Zeta assuming the underlying data X_1 and X_2 follow independent exponential distributions with different parameters.
#' @param z The coefficient of variation (CV) between two replicates
#' @param lam_1 The (positive) parameter for the first replicate's exponential distribution
#' @param lam_2 The (positive) parameter for the second replicate's exponential distribution
#' @return The marginal probability of observing Zeta <= z
#' @export
F_Z_asym <- function(z, lam_1, lam_2) {
  if(z < 0){
    return(0)
  }else if(z >= 0 && z <= sqrt(2)){
    return(
      lam_1*lam_2*(1-(g_helper(z))^2) / ( (lam_1+lam_2*g_helper(z)) * (lam_1*g_helper(z)+lam_2))
    )
  }else{
    return(1)
  }
}


#' Marginal CDF for Zeta Under Symmetric Exponential Assumption
#' 
#' @description This function is the marginal CDF of Zeta assuming the underlying data X_1 and X_2 follow independent and identically distributed exponential distributions.
#' @param z The coefficient of variation (CV) between two replicates
#' @param lam The (positive) parameter for the replicates's exponential distributions
#' @return The marginal probability of observing Zeta <= z
#' @export
F_Z_sym <- function(z, lam) {
  if(z < 0){
    return(0)
  }else if(z >= 0 && z <= sqrt(2)){
    return(
      (1-g_helper(z)) / (1+g_helper(z))
    )
  }else{
    return(1)
  }
}


#' Outlier Probability Under Asymmetric Exponential Assumption (Using Joint Method)
#' 
#' @description This function provides the outlier probability of a pair of replicates described by (Delta=d, Zeta=z) using the joint method assuming the underlying data X_1 and X_2 follow independent exponential distributions with different parameters.
#' @param d The absolute difference between two replicates
#' @param z The coefficient of variation (CV) between two replicates
#' @param lambda_1 The (positive) parameter for the first replicate's exponential distribution
#' @param lambda_2 The (positive) parameter for the second replicate's exponential distribution
#' @return The joint probability of observing Delta <= d and Z >= z (if d < 0) or Delta >= d and Z >= z (if d >= 0)
#' @export
q_exp_asym <- function(d, z, lambda_1, lambda_2) {
  if(z < 0 || z > sqrt(2)){
    return(1) #Must be some type of error to get this
  }else if(d < 0 && z >= 0 && z <= sqrt(2)){ #See picture of why this
    return(
      F_D_asym(d, lambda_1, lambda_2) - F_exp_asym(d, z, lambda_1, lambda_2)
    )
  }else if(d >= 0 && z >= 0 && z <= sqrt(2)){ #See picture of why this
    return(
      1 - F_D_asym(d, lambda_1, lambda_2) - F_Z_asym(z, lambda_1, lambda_2) +
        F_exp_asym(d, z, lambda_1, lambda_2)
    )
  }
}


#' Outlier Probability Under Symmetric Exponential Assumption (Using Joint Method)
#' 
#' @description This function provides the outlier probability of a pair of replicates described by (Delta=d, Zeta=z) using the joint method assuming the underlying data X_1 and X_2 follow independent and identically distributed exponential distributions.
#' @param d The absolute difference between two replicates
#' @param z The coefficient of variation (CV) between two replicates
#' @param lambda The (positive) parameter for the replicates's exponential distributions
#' @return The joint probability of observing Delta <= d and Z >= z (if d < 0) or Delta >= d and Z >= z (if d >= 0)
#' @export
q_exp_sym <- function(d, z, lambda) {
  if(z < 0 || z > sqrt(2)){
    return(1) #Must be some type of error to get this
  }else if(d < 0 && z >= 0 && z <= sqrt(2)){ #See picture of why this
    return(
      F_D_sym(d, lambda) - F_exp_sym(d, z, lambda)
    )
  }else if(d >= 0 && z >= 0 && z <= sqrt(2)){ #See picture of why this
    return(
      1 - F_D_sym(d, lambda) - F_Z_sym(z, lambda) + F_exp_sym(d, z, lambda)
    )
  }
}


#' Joint PDF Under GG Assumption
#' 
#' @description This function is the joint CDF of Delta and Zeta assuming the underlying data X_1 and X_2 follow independent (and possibly identically distributed) GG distributions. We use the parameterization used by the R package flexsurv, given by Kotz and Johnson (1970).
#' @param d The absolute difference between two replicates
#' @param z The coefficient of variation (CV) between two replicates
#' @param a_1 The parameter alpha > 0 in the GG distribution for X_1
#' @param b_1 The parameter beta > 0 in the GG distribution for X_1
#' @param c_1 The parameter c > 0 in the GG distribution for X_1
#' @param a_2 The parameter alpha > 0 in the GG distribution for X_2
#' @param b_2 The parameter beta > 0 in the GG distribution for X_2
#' @param c_2 The parameter c > 0 in the GG distribution for X_2
#' @return The joint PDF of observing Delta = d and Zeta = z, to be integrated to get useful probabilities
#' @export
f_gg <- function(d, z, a_1, b_1, c_1, a_2, b_2, c_2) {
  if(z < 0 || z > sqrt(2)){
    return(0)
  }else if(d >= 0 && z >= 0 && z <= sqrt(2)){
    return(
      (c_1/(b_1^(c_1*a_1)*gamma(a_1))) * ((sqrt(2)/z+1)*d/2)^(c_1*a_1-1) * exp(-((sqrt(2)/z+1)*d/(2*b_1))^c_1) *
        (c_2/(b_2^(c_2*a_2)*gamma(a_2))) * ((sqrt(2)/z-1)*d/2)^(c_2*a_2-1) * exp(-((sqrt(2)/z-1)*d/(2*b_2))^c_2) *
        sqrt(2)*d/(2*z^2)
    )
  }else if (d < 0 && z >= 0 && z <= sqrt(2)){
    return(
      (c_1/(b_1^(c_1*a_1)*gamma(a_1))) * ((-sqrt(2)/z+1)*d/2)^(c_1*a_1-1) * exp(-((-sqrt(2)/z+1)*d/(2*b_1))^c_1) *
        (c_2/(b_2^(c_2*a_2)*gamma(a_2))) * ((-sqrt(2)/z-1)*d/2)^(c_2*a_2-1) * exp(-((-sqrt(2)/z-1)*d/(2*b_2))^c_2) *
        (-sqrt(2)*d)/(2*z^2)
    )
  }
}


#' Outlier Probability Under GG Assumption (Using Joint Method)
#' 
#' @description This function provides the outlier probability of a pair of replicates described by (Delta=d, Zeta=z) using the joint method assuming the underlying data X_1 and X_2 follow independent (and possibly identically distributed) GG distributions. We use the parameterization used by the R package flexsurv, given by Kotz and Johnson (1970). We numerically integrate the PDF given by the function f_gg using the function adaptIntegrate in the package cubature.
#' @param d The absolute difference between two replicates
#' @param z The coefficient of variation (CV) between two replicates
#' @param alpha_1 The parameter alpha > 0 in the GG distribution for X_1
#' @param beta_1 The parameter beta > 0 in the GG distribution for X_1
#' @param c_1 The parameter c > 0 in the GG distribution for X_1
#' @param alpha_2 The parameter alpha > 0 in the GG distribution for X_2
#' @param beta_2 The parameter beta > 0 in the GG distribution for X_2
#' @param c_2 The parameter c > 0 in the GG distribution for X_2
#' @return The joint probability of observing Delta <= d and Z >= z (if d < 0) or Delta >= d and Z >= z (if d >= 0)
#' @import cubature
#' @export
q_gg <- function(d, z, alpha_1, beta_1, c_1, alpha_2, beta_2, c_2) {
  if(z < 0 || z > sqrt(2)){
    return(1) #Must be some type of error to get this
  }else if(d < 0 && z >= 0 && z <= sqrt(2)){
    return(
      adaptIntegrate((function(u) f_gg(d + 1 - 1/u[1], u[2], alpha_1, beta_1, c_1, alpha_2, beta_2, c_2)/(u[1]^2)),
                     c(0,z), c(1,sqrt(2)))$integral
    )
  }else if(d >= 0 && z >= 0 && z <= sqrt(2)){
    return(
      adaptIntegrate((function(u) f_gg(d + u[1]/(1-u[1]), u[2], alpha_1, beta_1, c_1, alpha_2, beta_2, c_2)/((1-u[1])^2)),
                     c(0,z), c(1,sqrt(2)))$integral
    )
  }
}


#' Outlier Probability for CV Under GG Assumption (Using Marginal Method)
#' 
#' @description This function provides the outlier probability of a pair of replicates described by Zeta = z using the marginal method assuming the underlying data X_1 and X_2 follow independent (and possibly identically distributed) GG distributions. We use the parameterization used by the R package flexsurv, given by Kotz and Johnson (1970). We numerically integrate the PDF given by the function f_gg using the function adaptIntegrate in the package cubature.
#' @param z The coefficient of variation (CV) between two replicates
#' @param alpha_1 The parameter alpha > 0 in the GG distribution for X_1
#' @param beta_1 The parameter beta > 0 in the GG distribution for X_1
#' @param c_1 The parameter c > 0 in the GG distribution for X_1
#' @param alpha_2 The parameter alpha > 0 in the GG distribution for X_2
#' @param beta_2 The parameter beta > 0 in the GG distribution for X_2
#' @param c_2 The parameter c > 0 in the GG distribution for X_2
#' @return The marginal probability of observing Zeta >= z
#' @import cubature
#' @export
q_gg_Z <- function(z, alpha_1, beta_1, c_1, alpha_2, beta_2, c_2){
  # Add together the contributions for D < 0 and D >= 0
  adaptIntegrate((function(u) f_gg(1 - 1/u[1], u[2], alpha_1, beta_1, c_1, alpha_2, beta_2, c_2)/(u[1]^2)),
                 c(0,z), c(1,sqrt(2)))$integral +
    adaptIntegrate((function(u) f_gg(u[1]/(1-u[1]), u[2], alpha_1, beta_1, c_1, alpha_2, beta_2, c_2)/((1-u[1])^2)),
                   c(0,z), c(1,sqrt(2)))$integral
}