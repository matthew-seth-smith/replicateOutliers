library(testthat)
library(flexsurv) #For rgengamma.orig

context("Test that outlier_DZ, q_exp_joint_DZ, and q_exp_marg_DZ work by recreating the Sim_GG data.frame.")

# Recreate the simulated GG data. This is the same procedure we used to make the GG data in the original paper
set.seed(1729) #Set the seed so we can replicate the random results, using Ramanujan's number for the seed
# I did not keep the seed from just generating the GG data, so I will need to re-simulate the gamma and
# Weibull data in order to have the seed for the GG data to be correct
n_total <- 10000 #Number of data points. Doing full simulation
# Gamma data:
df_gamma <- data.frame(
  section=c(rep(0, 0.2*n_total), rep(1, 0.39*n_total), rep(2, 0.39*n_total),
            rep(3, 0.01*n_total), rep(4, 0.01*n_total)),
  X_1=c(rexp(0.2*n_total, 1/5), rgamma(0.39*n_total, shape=120, rate=1),
        rgamma(0.39*n_total, shape=200, rate=1.5), rgamma(0.01*n_total, shape=100, rate=1),
        rgamma(0.01*n_total, shape=10, rate=1/2)),
  X_2=c(rexp(0.2*n_total, 1/6), rgamma(0.39*n_total, shape=200, rate= 1.5),
        rgamma(0.39*n_total, shape=120, rate=1), rgamma(0.01*n_total, shape=10, rate=1/2),
        rgamma(0.01*n_total, shape=100, rate=1))
)
# Weibull data:
df_wei <- data.frame(
  section=c(rep(0, 0.2*n_total), rep(1, 0.39*n_total), rep(2, 0.39*n_total),
            rep(3, 0.01*n_total), rep(4, 0.01*n_total)),
  X_1=c(rexp(0.2*n_total, 1/2), rweibull(0.39*n_total, shape=200, scale=150),
        rweibull(0.39*n_total, shape=10, scale=150), rweibull(0.01*n_total, shape=100, scale=170),
        rweibull(0.01*n_total, shape=30, scale=50)),
  X_2=c(rexp(0.2*n_total, 1/2), rweibull(0.39*n_total, shape=10, scale=150),
        rweibull(0.39*n_total, shape=200, scale=150), rweibull(0.01*n_total, shape=30, scale=50),
        rweibull(0.01*n_total, shape=100, scale=170))
)
# GG data:
df_gg <- data.frame(
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
df_gg$D <- df_gg$X_1 - df_gg$X_2
df_gg$Z <- sqrt(2) * abs(df_gg$X_1 - df_gg$X_2) / (df_gg$X_1 + df_gg$X_2)
df_gg$outlier_orig <- outlier_DZ(df_gg$D, df_gg$Z)
df_gg$q_exp_j <- q_exp_joint_DZ(df_gg$D, df_gg$Z)
df_gg$q_exp_m <- q_exp_marg_DZ(df_gg$D, df_gg$Z)

data(Sim_GG) #The data.frame contained in the package
test_that("The function outlier_DZ works as expected.",{
  expect_equivalent(df_gg$outlier_orig, Sim_GG$outlier_orig, tolerance=1e-5)
})
test_that("The function q_exp_joint_DZ works as expected.",{
  expect_equivalent(df_gg$q_exp_j, Sim_GG$q_exp_j, tolerance=1e-5)
})
test_that("The function q_exp_marg_DZ works as expected.",{
  expect_equivalent(df_gg$q_exp_m, Sim_GG$q_exp_m, tolerance=1e-5)
})