# Matthew Smith, 9/23/19
# In this file, we create and document all 5 of the outlier detection functions for this package


#' Original Asymmetric Laplace-Weibull Method for Outlier Detection
#' 
#' @description This function implements the original method for outlier detection among replicated data. It first fits the aboslute difference Delta between two replicates to an Asymmetric Laplace Distribution using MLE. Then among the points outside of some central band, we fit their values for the coefficient of variation Zeta to a Weibull (or lognormal) distribution (also using MLE). The points from this set with significantly large Zeta values are outliers.
#' @param D The absolute difference (Delta) between two vectors of (positive) replicated data: D = X_1 - X_2
#' @param Z The coefficient of variation (Zeta) between two vectors of (positive) replicated data: Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2)
#' @param type The distribution to which we fit Zeta. By default this value is "weibull", but we can also set it to "lognormal"
#' @param p_theta We use the (1-p_theta)*100\% two-sided confidence interval for theta in Delta = X_1 - X_2 + theta to determine if there is a significant translation of the absolute difference Delta. If this interval contains 0, then we set theta = 0. We set p_theta = 0.05 by default
#' @param p_kappa We use the (1-p_kappa)*100\% two-sided confidence interval for the asymmetry parameter kappa in the Asymmetric Laplace Distribution to which we fit Delta. If this interval for log(kappa) contains 0, then we set kappa = 0 and use a Symmetric Laplace Distribution for Delta. We set p_kappa = 0.05 by default
#' @param k The number of standard deviations about the center (mean) of the Asymmetric Laplace Distribution for Delta that we use to define the "central band." We set k = 1 by default
#' @param t We remove the top t*100\% of points when fitting the coefficient of variation Zeta values to a Weibull or lognormal distribution, so that these points do not cause problems when fitting the distribution. We set t = 0.01 by default
#' @param b We remove the bottom b*100\% of points when fitting the coefficient of variation Zeta values to a Weibull or lognormal distribution, so that these points do not cause problems when fitting the distribution. We set b = 0.01 by default
#' @param N_u After we fit the Zeta values for the correct subset of points to a Weibull or lognormal distribution, we return as outliers all the points whose Zeta values are greater than the level above which we would not expect to find N_u*100\% of the points. This procedure is how the underlying function getOutliersI in the package extremevalues works. We set N_u = 1 by default
#' @param N_l After we fit the Zeta values for the correct subset of points to a Weibull or lognormal distribution, we return as outliers all the points whose Zeta values are less than the level below which we would not expect to find N_l*100\% of the points. This procedure is how the underlying function getOutliersI in the package extremevalues works. We set N_l = 1 by default
#' @return A numerical vector of equal length to the input D and Z vectors, with values 0 for points in the central band, 1 for points in the "wings" with small Zeta, and 2 for outliers
#' @import VGAM
#' @import extremevalues
#' @examples
#' # Assume X_1 and X_2 are positive data vectors of the same length. These are the replicates
#' df <- data.frame(X_1=X_1, X_2=X_2)
#' df$D <- df$X_1 - df$X_2
#' df$Z <- sqrt(2) * abs(df$X_1 - df$X_2) / (df$X_1 + df$X_2)
#' df$outlier_status <- outlier_DZ(df$D, df$Z)
#' @export
outlier_DZ <- function(D, Z, type="weibull", p_theta=0.05,
                       p_kappa=0.05, k=1, t=0.01, b=0.01, N_u=0.01, N_l=0.01){
  if(length(D) != length(Z)){
    stop("The lengths of D and Z must be equal.")
  }
  
  # First fit D to Asymmetric Laplace
  fit_ala <- vglm(as.matrix(D) ~ 1, alaplace3, trace=TRUE, crit="c",
                  control = vglm.control(maxit=50))
  
  CI_theta <- Coef(fit_ala)[1] + summary(fit_ala)@coef3[1,2] * qnorm(1-p_theta/2) *
    c(-1,1) #95% confidence interval for theta using approximate normal for MLE
  CI_kappa <- log(Coef(fit_ala)[3]) + summary(fit_ala)@coef3[3,2] *
    qnorm(1-p_kappa/2) * c(-1,1)
  # 95% confidence interval for log(kappa) using approximate normal for MLE
  
  if(CI_theta[1] > 0 || CI_theta[2] < 0){ #If there IS a significant location parameter
    theta_est <- Coef(fit_ala)[1] #Make sure to subtract this later when plugging data into F_exp
  }else{
    theta_est <- 0 #This way we always have a value to subtract from D
  }
  # Estimate kappa, mu, and standard deviation
  if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If there IS significant asymmetry
    kappa_est <- Coef(fit_ala)[3]
    mu_est <- Coef(fit_ala)[2] * (1/kappa_est - kappa_est) / sqrt(2)
    sd_est <- sqrt(Coef(fit_ala)[2]^2 + mu_est^2)
  }else{
    kappa_est <- 1 #This way we always have a value of kappa, even if we won't use it
    mu_est <- 0 #Set equal to 0 to avoid numerical imprecision
    sd_est <- Coef(fit_ala)[2] #Set equal to this to avoid numerical imprecision
  }
  
  lap_range <- mu_est + theta_est + k * sd_est * c(-1,1) #Range of values to use for D
  
  df_temp <- data.frame(D=D, Z=Z, D_out=0, Z_out=0) #Initialize
  row.names(df_temp) <- 1:nrow(df_temp) #Explicitly label the rows
  df_temp$D_out <- as.numeric( #1 for TRUE and 0 for FALSE
    df_temp$D < lap_range[1] | df_temp$D > lap_range[2]
  )
  
  df_sub <- df_temp[df_temp$D_out==1,] #Subset of points with outlying D
  
  if(type=="weibull"){ #If using a Weibull distribution for Z
    ext_out <- getOutliersI(y=df_sub$Z, rho=c(N_u,N_l),
                            FLim=c(b,1-t), distribution="weibull")
  }else if(type=="lognormal"){
    ext_out <- getOutliersI(y=df_sub$Z, rho=c(N_u,N_l), 
                            FLim=c(b,1-t), distribution="lognormal")
  }else{
    stop("Distribution for Z must be \"weibull\" or \"lognormal\".")
  }
  
  df_sub$ext_out <- FALSE #Initialize
  df_sub$ext_out[ext_out$iLeft] <- TRUE #Small values of Z
  df_sub$ext_out[ext_out$iRight] <- TRUE #Large values of Z
  
  df_temp$Z_out[as.numeric(rownames(df_sub[df_sub$ext_out,]))] <- 1
  
  df_temp$D_out + df_temp$Z_out #Return combined measure of outlier status
}


#' Joint Exponential Method for Outlier Detection
#' 
#' @description This function implements the joint exponential method for outlier detection among replicated data. It first fits the aboslute difference Delta between two replicates to an Asymmetric Laplace Distribution using MLE. It then determines whether Delta's Laplace Distribution is Asymmetric or Symmetric and whether it has a significant displacement parameter. Using the exponential parameters derived from the Laplace parameters (see the paper in the citation), this function calculates the joint probability that Delta and Zeta take values more extreme than their observed values for each data point (d and z, respectively). That is, with Laplace displacement parameter theta (which can be zero), this function calculates P(D <= d, z <= Z <= sqrt(2)) if d < theta or P(D >= d, z <= Z <= sqrt(2)) else. We already know 0 <= Z <= sqrt(2).
#' @param D The absolute difference (Delta) between two vectors of (positive) replicated data: D = X_1 - X_2
#' @param Z The coefficient of variation (Zeta) between two vectors of (positive) replicated data: Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2)
#' @param p_theta We use the (1-p_theta)*100\% two-sided confidence interval for theta in Delta = X_1 - X_2 + theta to determine if there is a significant translation of the absolute difference Delta. If this interval contains 0, then we set theta = 0. We set p_theta = 0.05 by default
#' @param p_kappa We use the (1-p_kappa)*100\% two-sided confidence interval for the asymmetry parameter kappa in the Asymmetric Laplace Distribution to which we fit Delta. If this interval for log(kappa) contains 0, then we set kappa = 0 and use a Symmetric Laplace Distribution for Delta. We set p_kappa = 0.05 by default
#' @return A numerical vector of equal length to the input D and Z vectors, with the outlier probability q = P((D <= d, z <= Z <= sqrt(2)) | d < theta OR (D >= d, z <= Z <= sqrt(2)) | d >= theta)
#' @import VGAM
#' @examples
#' # Assume X_1 and X_2 are positive data vectors of the same length. These are the replicates
#' df <- data.frame(X_1=X_1, X_2=X_2)
#' df$D <- df$X_1 - df$X_2
#' df$Z <- sqrt(2) * abs(df$X_1 - df$X_2) / (df$X_1 + df$X_2)
#' df$q_exp_j <- q_exp_joint_DZ(df$D, df$Z)
#' @export
q_exp_joint_DZ <- function(D, Z, p_theta=0.05, p_kappa=0.05){
  if(length(D) != length(Z)){
    stop("The lengths of D and Z must be equal.")
  }
  
  # First fit D to Asymmetric Laplace
  fit_ala <- vglm(as.matrix(D) ~ 1, alaplace3, trace=TRUE, crit="c",
                  control = vglm.control(maxit=50))
  
  CI_theta <- Coef(fit_ala)[1] + summary(fit_ala)@coef3[1,2] * qnorm(1-p_theta/2) *
    c(-1,1) #95% confidence interval for theta using approximate normal for MLE
  CI_kappa <- log(Coef(fit_ala)[3]) + summary(fit_ala)@coef3[3,2] *
    qnorm(1-p_kappa/2) * c(-1,1)
  # 95% confidence interval for log(kappa) using approximate normal for MLE
  
  if(CI_theta[1] > 0 || CI_theta[2] < 0){ #If there IS a significant location parameter
    theta_est <- Coef(fit_ala)[1] #Make sure to subtract this later when plugging data into F_exp
  }else{
    theta_est <- 0 #This way we always have a value to subtract from D
  }
  
  if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If asymmetric. Use q_exp_asym
    lambda_exp_1est <- sqrt(2)*(Coef(fit_ala)[3])/(Coef(fit_ala)[2])
    lambda_exp_2est <- sqrt(2)/(Coef(fit_ala)[2]*Coef(fit_ala)[3])
  }else{ #If symmetric. Use q_exp_sym
    lambda_exp_est <- sqrt(2)/(Coef(fit_ala)[2])
  }
  
  # Probability calculation
  if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If we determine the data is asymmetric
    return(
      mapply(function(d, z)
      {q_exp_asym(d, z, lambda_exp_1est, lambda_exp_2est)}, D - theta_est, Z
      )
    )
  }else{ #If the data is symmetric
    return(
      mapply(function(d, z)
      {q_exp_sym(d, z, lambda_exp_est)}, D - theta_est, Z
      )
    )
  }
}


#' Marginal Exponential Method for Outlier Detection
#' 
#' @description This function implements the marginal exponential method for outlier detection among replicated data. It first fits the aboslute difference Delta between two replicates to an Asymmetric Laplace Distribution using MLE. It then determines whether Delta's Laplace Distribution is Asymmetric or Symmetric and whether it has a significant displacement parameter. Then among the points outside of some central band, we use the exponential parameters fitted to the Laplace Distribution for Delta (see the paper in the citation) to determine the marginal probability that Z will take a value greater than its observed value. We assign the probability 1 to points in the middle band.
#' @param D The absolute difference (Delta) between two vectors of (positive) replicated data: D = X_1 - X_2
#' @param Z The coefficient of variation (Zeta) between two vectors of (positive) replicated data: Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2)
#' @param p_theta We use the (1-p_theta)*100\% two-sided confidence interval for theta in Delta = X_1 - X_2 + theta to determine if there is a significant translation of the absolute difference Delta. If this interval contains 0, then we set theta = 0. We set p_theta = 0.05 by default
#' @param p_kappa We use the (1-p_kappa)*100\% two-sided confidence interval for the asymmetry parameter kappa in the Asymmetric Laplace Distribution to which we fit Delta. If this interval for log(kappa) contains 0, then we set kappa = 0 and use a Symmetric Laplace Distribution for Delta. We set p_kappa = 0.05 by default
#' @param k The number of standard deviations about the center (mean) of the Asymmetric Laplace Distribution for Delta that we use to define the "central band." We set k = 1 by default
#' @return A numerical vector of equal length to the input D and Z vectors, with the outlier probability q = P(z <= Z <= sqrt(2)) if (d,z) is not in the middle band and the assigned value 1 if it is in the middle band
#' @import VGAM
#' @examples
#' # Assume X_1 and X_2 are positive data vectors of the same length. These are the replicates
#' df <- data.frame(X_1=X_1, X_2=X_2)
#' df$D <- df$X_1 - df$X_2
#' df$Z <- sqrt(2) * abs(df$X_1 - df$X_2) / (df$X_1 + df$X_2)
#' df$q_exp_m <- q_exp_marg_DZ(df$D, df$Z)
#' @export
q_exp_marg_DZ <- function(D, Z, p_theta=0.05, p_kappa=0.05, k=1){
  if(length(D) != length(Z)){
    stop("The lengths of D and Z must be equal.")
  }
  
  # Assign the value 1 to the middle band in the marginal method
  # First fit D to Asymmetric Laplace
  fit_ala <- vglm(as.matrix(D) ~ 1, alaplace3, trace=TRUE, crit="c",
                  control = vglm.control(maxit=50))
  
  CI_theta <- Coef(fit_ala)[1] + summary(fit_ala)@coef3[1,2] * qnorm(1-p_theta/2) *
    c(-1,1) #95% confidence interval for theta using approximate normal for MLE
  CI_kappa <- log(Coef(fit_ala)[3]) + summary(fit_ala)@coef3[3,2] *
    qnorm(1-p_kappa/2) * c(-1,1)
  # 95% confidence interval for log(kappa) using approximate normal for MLE
  
  if(CI_theta[1] > 0 || CI_theta[2] < 0){ #If there IS a significant location parameter
    theta_est <- Coef(fit_ala)[1] #Make sure to subtract this later when plugging data into F_exp
  }else{
    theta_est <- 0 #This way we always have a value to subtract from D
  }
  # Estimate kappa, mu, and standard deviation
  if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If there IS significant asymmetry
    kappa_est <- Coef(fit_ala)[3]
    mu_est <- Coef(fit_ala)[2] * (1/kappa_est - kappa_est) / sqrt(2)
    sd_est <- sqrt(Coef(fit_ala)[2]^2 + mu_est^2)
    lambda_exp_1est <- sqrt(2)*(Coef(fit_ala)[3])/(Coef(fit_ala)[2])
    lambda_exp_2est <- sqrt(2)/(Coef(fit_ala)[2]*Coef(fit_ala)[3])
  }else{
    kappa_est <- 1 #This way we always have a value of kappa, even if we won't use it
    mu_est <- 0 #Set equal to 0 to avoid numerical imprecision
    sd_est <- Coef(fit_ala)[2] #Set equal to this to avoid numerical imprecision
    lambda_exp_est <- sqrt(2)/(Coef(fit_ala)[2])
  }
  
  lap_range <- mu_est + theta_est + k * sd_est * c(-1,1) #Range of values to use for D
  
  df_temp <- data.frame(D=D, Z=Z, D_out=0, Z_q=1) #Initialize
  row.names(df_temp) <- 1:nrow(df_temp) #Explicitly label the rows
  df_temp$D_out <- as.numeric( #1 for TRUE and 0 for FALSE
    df_temp$D < lap_range[1] | df_temp$D > lap_range[2]
  )
  
  if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If there IS significant asymmetry
    df_temp$Z_q[df_temp$D_out==1] <-
      sapply(df_temp$Z[df_temp$D_out==1], function(z){
        1 - F_Z_asym(z, lambda_exp_1est, lambda_exp_2est)
      })
  }else{
    df_temp$Z_q[df_temp$D_out==1] <-
      sapply(df_temp$Z[df_temp$D_out==1], function(z){
        1 - F_Z_sym(z, lambda_exp_est)
      })
  }
  
  df_temp$Z_q #Return the marginal q value
}


#' Joint Generalized Gamma Method for Outlier Detection
#' 
#' @description This function implements the joint generalized gamma method for outlier detection among replicated data. It first fits each replicate (X_1 and X_2) to generalized gamma distribution (using the parameterization in the R package flexsurv, given by Kotz and Johnson (1970)) using MLE. We then use those parameters and the joint probability density function derived in the paper in the citation for this package to calculate the joint probability that D and Z take values more extreme than their observed values for each data point (d and z, respectively). That is, we use numerical integration (using the function adaptIntegrate in the package cubature) the joint PDF to find P(D <= d, z <= Z <= sqrt(2)) if d < 0 or P(D >= d, z <= Z <= sqrt(2)) else. We already know 0 <= Z <= sqrt(2).
#' @param X_1 The first (independent) replicate of the data. A vector of positive real numbers
#' @param X_2 The second (independent) replicate of the data. A vector of positive real numbers
#' @param n_cores This function works by numerically integrating the joint PDF for each data point. To speed up this process, we run this process in parallel (using the package parallel), which requires specifying the number of cores (n_cores) on the computer to use. By default, we use all but one core on the machine (with the remaining one free for other functions).
#' @return A numerical vector of equal length to the input X_1 and X_2 vectors. Using D = X_1 - X_2, Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2), and (d,z) for each (X_1,X_2) data point, we get the outlier probability q = P((D <= d, z <= Z <= sqrt(2)) | d < theta OR (D >= d, z <= Z <= sqrt(2)) | d >= theta) for each (d,z)
#' @import flexsurv
#' @import cubature
#' @import parallel
#' @examples
#' # Assume X_1 and X_2 are positive data vectors of the same length. These are the replicates
#' df <- data.frame(X_1=X_1, X_2=X_2)
#' # The function q_gg_joint_DZ calculates D and Z for us
#' df$q_gg_j <- q_gg_joint_DZ(df$X_1, df$X_2)
#' @export
q_gg_joint_DZ <- function(X_1, X_2, n_cores=detectCores()-1){ #Joint method
  # This function starts with X_1 and X_2, not D and Z. We call it DZ because
  # we use those variables for the probability calculations
  
  if(length(X_1) != length(X_2)){
    stop("The lengths of X_1 and X_2 must be equal.")
  }
  
  # MLE
  fit_gg_1 <- flexsurvreg(Surv(X_1) ~ 1, dist="gengamma.orig")
  coef_1 <- exp(fit_gg_1$coefficients) #Need to exponentiate
  alpha_gg_1est <- coef_1[3]
  beta_1est <- coef_1[2]
  c_gg_1est <- coef_1[1]
  fit_gg_2 <- flexsurvreg(Surv(X_2) ~ 1, dist="gengamma.orig")
  coef_2 <- exp(fit_gg_2$coefficients) #Need to exponentiate
  alpha_gg_2est <- coef_2[3]
  beta_2est <- coef_2[2]
  c_gg_2est <- coef_2[1]
  # Test to see if actually Weibull (alpha=1) or actually gamma (c=1):
  # Can't figure out how to extract standard error from fit_gg object for confidence interval. Figure out later.
  
  D <- X_1 - X_2 #Find the difference
  Z <- sqrt(2)*abs(X_1-X_2)/(X_1+X_2) #Find the coefficient of variation
  
  # Probability calculation
  cl <- makeCluster(n_cores, type="FORK")
  return(
    unlist(parLapply(
      cl, 1:length(D), function(i){
        q_gg(D[i], Z[i], alpha_gg_1est, beta_1est, c_gg_1est, alpha_gg_2est, beta_2est, c_gg_2est)
      }
    ), use.names=FALSE)
  )
  stopCluster(cl)
}


#' Marginal Generalized Gamma Method for Outlier Detection
#' 
#' @description This function implements the marginal generalized gamma method for outlier detection among replicated data. It first fits each replicate (X_1 and X_2) to generalized gamma distributions (using the parameterization in the R package flexsurv, given by Kotz and Johnson (1970)) using MLE. It also fits the aboslute difference Delta (D = X_1 - X_2) between the two replicates to an Asymmetric Laplace Distribution using MLE. It then determines whether Delta's Laplace Distribution is Asymmetric or Symmetric and whether it has a significant displacement parameter. Then among the points outside of some central band (defined using the Laplace parameters fitted to Delta), we use the generalized gamma parameters fitted to the entire X_1 and X_2 vectors (see the paper in the citation) to determine the marginal probability that Z will take a value greater than its observed value. We use numerical integration (specifically the function adaptIntegrate in the package cubature) to integrate the marginal PDF for Z to get this probability. We assign the probability 1 to points in the middle band.
#' @param X_1 The first (independent) replicate of the data. A vector of positive real numbers
#' @param X_2 The second (independent) replicate of the data. A vector of positive real numbers
#' @param p_theta We use the (1-p_theta)*100\% two-sided confidence interval for theta in Delta = X_1 - X_2 + theta to determine if there is a significant translation of the absolute difference Delta. If this interval contains 0, then we set theta = 0. We set p_theta = 0.05 by default
#' @param p_kappa We use the (1-p_kappa)*100\% two-sided confidence interval for the asymmetry parameter kappa in the Asymmetric Laplace Distribution to which we fit Delta. If this interval for log(kappa) contains 0, then we set kappa = 0 and use a Symmetric Laplace Distribution for Delta. We set p_kappa = 0.05 by default
#' @param k The number of standard deviations about the center (mean) of the Asymmetric Laplace Distribution for Delta that we use to define the "central band." We set k = 1 by default
#' @param n_cores This function works by numerically integrating the joint PDF for each data point. To speed up this process, we run this process in parallel (using the package parallel), which requires specifying the number of cores (n_cores) on the computer to use. By default, we use all but one core on the machine (with the remaining one free for other functions).
#' @return A numerical vector of equal length to the input X_1 and X_2 vectors. Using D = X_1 - X_2, Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2), and (d,z) for each (X_1,X_2) data point, we get the marginal probability q = P(z <= Z <= sqrt(2)) if (d,z) is not in the middle band and the assigned value 1 if it is in the middle band
#' @import flexsurv
#' @import VGAM
#' @import cubature
#' @import parallel
#' @examples
#' # Assume X_1 and X_2 are positive data vectors of the same length. These are the replicates
#' df <- data.frame(X_1=X_1, X_2=X_2)
#' # The function q_gg_marg_DZ calculates D and Z for us
#' df$q_gg_m <- q_gg_marg_DZ(df$X_1, df$X_2)
#' @export
q_gg_marg_DZ <- function(X_1, X_2, p_theta=0.05, p_kappa=0.05, k=1, n_cores=detectCores()-1){
  # This function starts with X_1 and X_2, not D and Z. We call it DZ because
  # we use those variables for the probability calculations
  
  D <- X_1 - X_2 #Find the difference
  Z <- sqrt(2)*abs(X_1-X_2)/(X_1+X_2) #Find the coefficient of variation
  
  # MLE for X_1 and X_2
  fit_gg_1 <- flexsurvreg(Surv(X_1) ~ 1, dist="gengamma.orig")
  coef_1 <- exp(fit_gg_1$coefficients) #Need to exponentiate
  alpha_gg_1est <- coef_1[3]
  beta_1est <- coef_1[2]
  c_gg_1est <- coef_1[1]
  fit_gg_2 <- flexsurvreg(Surv(X_2) ~ 1, dist="gengamma.orig")
  coef_2 <- exp(fit_gg_2$coefficients) #Need to exponentiate
  alpha_gg_2est <- coef_2[3]
  beta_2est <- coef_2[2]
  c_gg_2est <- coef_2[1]
  # Test to see if actually Weibull (alpha=1) or actually gamma (c=1):
  # Can't figure out how to extract standard error from fit_gg object for confidence interval. Figure out later.
  
  # MLE for D
  fit_ala <- vglm(as.matrix(D) ~ 1, alaplace3, trace=TRUE, crit="c", control = vglm.control(maxit=50))
  
  CI_theta <- Coef(fit_ala)[1] + summary(fit_ala)@coef3[1,2] * qnorm(1-p_theta/2) *
    c(-1,1) #95% confidence interval for theta using approximate normal for MLE
  CI_kappa <- log(Coef(fit_ala)[3]) + summary(fit_ala)@coef3[3,2] *
    qnorm(1-p_kappa/2) * c(-1,1)
  # 95% confidence interval for log(kappa) using approximate normal for MLE
  
  if(CI_theta[1] > 0 || CI_theta[2] < 0){ #If there IS a significant location parameter
    theta_est <- Coef(fit_ala)[1]
  }else{
    theta_est <- 0 #This way we always have a value to subtract from D
  }
  
  # Estimate kappa, mu, and standard deviation
  if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If there IS significant asymmetry
    kappa_est <- Coef(fit_ala)[3]
    mu_est <- Coef(fit_ala)[2] * (1/kappa_est - kappa_est) / sqrt(2)
    sd_est <- sqrt(Coef(fit_ala)[2]^2 + mu_est^2)
  }else{
    kappa_est <- 1 #This way we always have a value of kappa, even if we won't use it
    mu_est <- 0 #Set equal to 0 to avoid numerical imprecision
    sd_est <- Coef(fit_ala)[2] #Set equal to this to avoid numerical imprecision
  }
  
  lap_range <- mu_est + theta_est + k * sd_est * c(-1,1) #Range of values to use for D
  
  df_temp <- data.frame(D=D, Z=Z, D_out=0, Z_q=1) #Initialize
  row.names(df_temp) <- 1:nrow(df_temp) #Explicitly label the rows
  df_temp$D_out <- as.numeric( #1 for TRUE and 0 for FALSE
    df_temp$D < lap_range[1] | df_temp$D > lap_range[2]
  )
  
  # Probability calculation
  cl <- makeCluster(n_cores, type="FORK")
  df_temp$Z_q[df_temp$D_out==1] <- unlist(
    parLapply(
      cl, which(df_temp$D_out==1), function(i){
        q_gg_Z(Z[i], alpha_gg_1est, beta_1est, c_gg_1est, alpha_gg_2est, beta_2est, c_gg_2est)
      }
    ), use.names=FALSE
  )
  stopCluster(cl)
  
  df_temp$Z_q #Return the marginal q value
}