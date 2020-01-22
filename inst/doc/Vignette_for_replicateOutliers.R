## ---- message=FALSE, warning=FALSE---------------------------------------
rm(list=ls())
set.seed(1729)
# Set the seed so we can replicate the random results, using Ramanujan's number for the seed
# I removed some of the steps I used in making the original Sim_GG data.frame, so this seed will not
# exactly match the one used to make Sim_GG

# library(devtools) #For document
# document()
# library(replicateOutliers) #To get the functions for these methods
library(flexsurv) #For GG functions
library(ggplot2) #For ggplot plotting functions

n_total <- 10000 #Number of data points. Doing full simulation
df <- data.frame(
  section=c(rep(0, 0.2*n_total), rep(1, 0.39*n_total), rep(2, 0.39*n_total),
            rep(3, 0.01*n_total), rep(4, 0.01*n_total)),
  X_1=c(rexp(0.2*n_total, 1/5), rgengamma.orig(0.39*n_total, shape=1.1, scale=1.6, k=120),
        rgengamma.orig(0.39*n_total, shape=0.95, scale=2/3, k=200),
        rgengamma.orig(0.01*n_total, shape=0.898, scale=2/3, k=100),
        rgengamma.orig(0.01*n_total, shape=0.89, scale=1.6, k=125)),
  X_2=c(rexp(0.2*n_total, 1/5), rgengamma.orig(0.39*n_total, shape=0.95, scale=2/3, k=200),
        rgengamma.orig(0.39*n_total, shape=1.1, scale=1.6, k=120),
        rgengamma.orig(0.01*n_total, shape=0.89, scale=1.6, k=125),
        rgengamma.orig(0.01*n_total, shape=0.898, scale=2/3, k=100))
)
df$D <- df$X_1 - df$X_2
df$Z <- sqrt(2) * abs(df$X_1 - df$X_2) / (df$X_1 + df$X_2)

head(df)

ggplot(df, aes(x=D, y=Z, color=section)) + geom_point() +
  ggtitle("Volcano Plot Using Generalized Gamma Distribution\n(New Simulated Data)")

## ------------------------------------------------------------------------
data(Sim_GG)
head(Sim_GG)
ggplot(df, aes(x=D, y=Z, color=section)) + geom_point() +
  ggtitle("Volcano Plot Using Generalized Gamma Distribution\n(Sim_GG Data Frame)")

## ------------------------------------------------------------------------
library(viridis) #For scale_colour_viridis in ggplot
df$q_exp_joint <- q_exp_joint_DZ(df$D, df$Z)

ggplot(df, aes(x=D, y=Z, color=q_exp_joint)) + scale_colour_viridis() +
  geom_point() + labs(color="q-Value") + ggtitle("Joint Exponential Method")

## ------------------------------------------------------------------------
library(gridExtra) #For grid.arrange
ej_1 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_joint<0.05))) +
  geom_point() + labs(color=expression(q<0.05)) + ggtitle("Joint Exponential")
ej_2 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_joint<0.01))) +
  geom_point() + labs(color=expression(q<0.01)) + ggtitle("Joint Exponential")
ej_3 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_joint<0.005))) +
  geom_point() + labs(color=expression(q<0.005)) + ggtitle("Joint Exponential")
ej_4 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_joint<0.001))) +
  geom_point() + labs(color=expression(q<0.001)) + ggtitle("Joint Exponential")
ej_5 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_joint<0.0005))) +
  geom_point() + labs(color=expression(q<5%*%10^-4)) + ggtitle("Joint Exponential")
ej_6 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_joint<0.0001))) +
  geom_point() + labs(color=expression(q<1%*%10^-4)) + ggtitle("Joint Exponential")
grid.arrange(ej_1, ej_2, ej_3, ej_4, ej_5, ej_6, nrow=3, ncol=2)


sum(df$q_exp_joint < 0.05) / n_total #0.0415
sum(df$q_exp_joint < 0.01) / n_total #0.0232
sum(df$q_exp_joint < 0.005) / n_total #0.0182
sum(df$q_exp_joint < 0.001) / n_total #0.0035
sum(df$q_exp_joint < 5e-4) / n_total #7e-04
sum(df$q_exp_joint < 1e-4) / n_total #0

## ------------------------------------------------------------------------
df$q_exp_marg <- q_exp_marg_DZ(df$D, df$Z)

ggplot(df, aes(x=D, y=Z, color=q_exp_marg)) + scale_colour_viridis() +
  geom_point() + labs(color="q-Value") + ggtitle("Marginal Exponential Method")

mj_1 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_marg<0.7))) +
  geom_point() + labs(color=expression(q<0.7)) + ggtitle("Marginal Exponential")
mj_2 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_marg<0.6))) +
  geom_point() + labs(color=expression(q<0.6)) + ggtitle("Marginal Exponential")
mj_3 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_marg<0.5))) +
  geom_point() + labs(color=expression(q<0.5)) + ggtitle("Marginal Exponential")
mj_4 <- ggplot(df, aes(x=D, y=Z, color=(q_exp_marg<0.4))) +
  geom_point() + labs(color=expression(q<0.4)) + ggtitle("Marginal Exponential")
grid.arrange(mj_1, mj_2, mj_3, mj_4, nrow=2, ncol=2)

sum(df$q_exp_marg < 0.7) / n_total #0.0263
sum(df$q_exp_marg < 0.6) / n_total #0.0199
sum(df$q_exp_marg < 0.5) / n_total #0.0131
sum(df$q_exp_marg < 0.4) / n_total #0.0012

## ---- eval=FALSE---------------------------------------------------------
#  #' Original Asymmetric Laplace-Weibull Method for Outlier Detection
#  #'
#  #' @description This function implements the original method for outlier detection among replicated data. It first fits the aboslute difference Delta between two replicates to an Asymmetric Laplace Distribution using MLE. Then among the points outside of some central band, we fit their values for the coefficient of variation Zeta to a Weibull (or lognormal) distribution (also using MLE). The points from this set with significantly large Zeta values are outliers.
#  #' @param D The absolute difference (Delta) between two vectors of (positive) replicated data: D = X_1 - X_2
#  #' @param Z The coefficient of variation (Zeta) between two vectors of (positive) replicated data: Z = sqrt(2) * abs(X_1 - X_2) / (X_1 + X_2)
#  #' @param type The distribution to which we fit Zeta. By default this value is "weibull", but we can also set it to "lognormal"
#  #' @param p_theta We use the (1-p_theta)*100\% two-sided confidence interval for theta in Delta = X_1 - X_2 + theta to determine if there is a significant translation of the absolute difference Delta. If this interval contains 0, then we set theta = 0. We set p_theta = 0.05 by default
#  #' @param p_kappa We use the (1-p_kappa)*100\% two-sided confidence interval for the asymmetry parameter kappa in the Asymmetric Laplace Distribution to which we fit Delta. If this interval for log(kappa) contains 0, then we set kappa = 0 and use a Symmetric Laplace Distribution for Delta. We set p_kappa = 0.05 by default
#  #' @param k The number of standard deviations about the center (mean) of the Asymmetric Laplace Distribution for Delta that we use to define the "central band." We set k = 1 by default
#  #' @param t We remove the top t*100\% of points when fitting the coefficient of variation Zeta values to a Weibull or lognormal distribution, so that these points do not cause problems when fitting the distribution. We set t = 0.01 by default
#  #' @param b We remove the bottom b*100\% of points when fitting the coefficient of variation Zeta values to a Weibull or lognormal distribution, so that these points do not cause problems when fitting the distribution. We set b = 0.01 by default
#  #' @param N_u After we fit the Zeta values for the correct subset of points to a Weibull or lognormal distribution, we return as outliers all the points whose Zeta values are greater than the level above which we would not expect to find N_u*100\% of the points. This procedure is how the underlying function getOutliersI in the package extremevalues works. We set N_u = 1 by default
#  #' @param N_l After we fit the Zeta values for the correct subset of points to a Weibull or lognormal distribution, we return as outliers all the points whose Zeta values are less than the level below which we would not expect to find N_l*100\% of the points. This procedure is how the underlying function getOutliersI in the package extremevalues works. We set N_l = 1 by default
#  #' @return A numerical vector of equal length to the input D and Z vectors, with values 0 for points in the central band, 1 for points in the "wings" with small Zeta, and 2 for outliers
#  #' @import VGAM
#  #' @import extremevalues
#  #' @examples
#  #' # Assume X_1 and X_2 are positive data vectors of the same length. These are the replicates
#  #' df <- data.frame(X_1=Sim_GG$X_1, X_2=Sim_GG$X_2)
#  #' df$D <- df$X_1 - df$X_2
#  #' df$Z <- sqrt(2) * abs(df$X_1 - df$X_2) / (df$X_1 + df$X_2)
#  #' df$outlier_status <- outlier_DZ(df$D, df$Z)
#  #' @export
#  outlier_DZ <- function(D, Z, type="weibull", p_theta=0.05,
#                         p_kappa=0.05, k=1, t=0.01, b=0.01, N_u=0.01, N_l=0.01){
#    if(length(D) != length(Z)){
#      stop("The lengths of D and Z must be equal.")
#    }
#  
#    # First fit D to Asymmetric Laplace
#    fit_ala <- vglm(as.matrix(D) ~ 1, alaplace3, trace=TRUE, crit="c",
#                    control = vglm.control(maxit=50))
#  
#    CI_theta <- Coef(fit_ala)[1] + summary(fit_ala)@coef3[1,2] * qnorm(1-p_theta/2) *
#      c(-1,1) #95% confidence interval for theta using approximate normal for MLE
#    CI_kappa <- log(Coef(fit_ala)[3]) + summary(fit_ala)@coef3[3,2] *
#      qnorm(1-p_kappa/2) * c(-1,1)
#    # 95% confidence interval for log(kappa) using approximate normal for MLE
#  
#    if(CI_theta[1] > 0 || CI_theta[2] < 0){ #If there IS a significant location parameter
#      theta_est <- Coef(fit_ala)[1] #Make sure to subtract this later when plugging data into F_exp
#    }else{
#      theta_est <- 0 #This way we always have a value to subtract from D
#    }
#    # Estimate kappa, mu, and standard deviation
#    if(CI_kappa[1] > 0 || CI_kappa[2] < 0){ #If there IS significant asymmetry
#      kappa_est <- Coef(fit_ala)[3]
#      mu_est <- Coef(fit_ala)[2] * (1/kappa_est - kappa_est) / sqrt(2)
#      sd_est <- sqrt(Coef(fit_ala)[2]^2 + mu_est^2)
#    }else{
#      kappa_est <- 1 #This way we always have a value of kappa, even if we won't use it
#      mu_est <- 0 #Set equal to 0 to avoid numerical imprecision
#      sd_est <- Coef(fit_ala)[2] #Set equal to this to avoid numerical imprecision
#    }
#  
#    lap_range <- mu_est + theta_est + k * sd_est * c(-1,1) #Range of values to use for D
#  
#    df_temp <- data.frame(D=D, Z=Z, D_out=0, Z_out=0) #Initialize
#    row.names(df_temp) <- 1:nrow(df_temp) #Explicitly label the rows
#    df_temp$D_out <- as.numeric( #1 for TRUE and 0 for FALSE
#      df_temp$D < lap_range[1] | df_temp$D > lap_range[2]
#    )
#  
#    df_sub <- df_temp[df_temp$D_out==1,] #Subset of points with outlying D
#  
#    if(type=="weibull"){ #If using a Weibull distribution for Z
#      ext_out <- getOutliersI(y=df_sub$Z, rho=c(N_u,N_l),
#                              FLim=c(b,1-t), distribution="weibull")
#    }else if(type=="lognormal"){
#      ext_out <- getOutliersI(y=df_sub$Z, rho=c(N_u,N_l),
#                              FLim=c(b,1-t), distribution="lognormal")
#    }else{
#      stop("Distribution for Z must be \"weibull\" or \"lognormal\".")
#    }
#  
#    df_sub$ext_out <- FALSE #Initialize
#    df_sub$ext_out[ext_out$iLeft] <- TRUE #Small values of Z
#    df_sub$ext_out[ext_out$iRight] <- TRUE #Large values of Z
#  
#    df_temp$Z_out[as.numeric(rownames(df_sub[df_sub$ext_out,]))] <- 1
#  
#    df_temp$D_out + df_temp$Z_out #Return combined measure of outlier status
#  }

## ------------------------------------------------------------------------
citation("devtools")

citation("flexsurv")

citation("ggplot2")

citation("VGAM")

citation("viridis")

citation("gridExtra")

library(cubature)
citation("cubature")

library(parallel)
citation("parallel")

