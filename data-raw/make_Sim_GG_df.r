# Matthew Smith, 10/6/19
# Use this file to import the simulated GG data set
setwd("~/Documents/FCCC_Research/Ross_Lab/New_Stuff/replicateOutliers/data-raw")
rm(list=ls())
set.seed(1729) #Set the seed so we can replicate any random results, using
# Ramanujan's number for the seed


library(devtools)
document()
library(replicateOutliers)


df <- read.csv("Simulated_GG.csv", header=TRUE)
df$Z <- df$CV #Rename
df$q_gg_j <- df$gg_j #Rename
df$q_gg_m <- df$gg_m #Rename
# This already has the GG methods, so we just need to do the original
# Asymmetric Laplace/Weibull method and the exponential methods
df$outlier_orig <- outlier_DZ(df$D, df$Z)
df$q_exp_j <- q_exp_joint_DZ(df$D, df$Z)
df$q_exp_m <- q_exp_marg_DZ(df$D, df$Z)


# Keep the "section" variable in case we want to use it
Sim_GG <- df[,c("section", "X_1", "X_2", "D", "Z", "outlier_orig",
                "q_exp_j", "q_exp_m", "q_gg_j", "q_gg_m")]

save(Sim_GG, file = "../data/Sim_GG.rda")
