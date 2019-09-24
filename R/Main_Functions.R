# Matthew Smith, 9/23/19
# In this file, we create and document all 5 of the outlier detection functions for this package


#' Original Asymmetric Laplace-Weibull Method for Outlier Detection
#' 
#' @description This function implements the original method for outlier detection. It first fits the aboslute difference Delta between two replicates to an Asymmetric Laplace Distribution using MLE. Then among the points outside of some central band, we fit their values for the coefficient of variation Zeta to a Weibull (or lognormal) distribution (also using MLE). The points from this set with significantly large Zeta values are outliers.
#' @param D The absolute difference (Delta) between two replicates: D = X_1 - X_2
#' @param Z The coefficient of variation (Zeta) between two replicates: Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2)
#' @param type The distribution to which we fit Zeta. By default this value is "weibull", but we can also set it to "lognormal"
#' @param p_theta We use the (1-p_theta)*100\% two-sided confidence interval for theta in Delta = X_1 - X_2 + theta to determine if there is a significant translation of the absolute difference Delta. If this interval contains 0, then we set theta = 0. We set p_theta = 0.05 by default
#' @param p_kappa We use the (1-p_kappa)*100\% two-sided confidence interval for the asymmetry parameter kappa in the Asymmetric Laplace Distribution to which we fit Delta. If this interval for log(kappa) contains 0, then we set kappa = 0 and use a Symmetric Laplace Distribution for Delta. We set p_kappa = 0.05 by default
#' @param k The number of standard deviations about the center of the Asymmetric Laplace Distribution for Delta that we use to define the "central band." We set k = 1 by default
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
    stop("Distribution for CV must be \"weibull\" or \"lognormal\".")
  }
  
  df_sub$ext_out <- FALSE #Initialize
  df_sub$ext_out[ext_out$iLeft] <- TRUE #Small values of Z
  df_sub$ext_out[ext_out$iRight] <- TRUE #Large values of Z
  
  df_temp$Z_out[as.numeric(rownames(df_sub[df_sub$ext_out,]))] <- 1
  
  df_temp$D_out + df_temp$Z_out #Return combined measure of outlier status
}