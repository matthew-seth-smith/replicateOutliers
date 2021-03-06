% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helper_Functions.R
\name{F_Z_asym}
\alias{F_Z_asym}
\title{Marginal CDF for Zeta Under Asymmetric Exponential Assumption}
\usage{
F_Z_asym(z, lam_1, lam_2)
}
\arguments{
\item{z}{The coefficient of variation (CV) between two replicates}

\item{lam_1}{The (positive) parameter for the first replicate's exponential distribution}

\item{lam_2}{The (positive) parameter for the second replicate's exponential distribution}
}
\value{
The marginal probability of observing Zeta <= z
}
\description{
This function is the marginal CDF of Zeta assuming the underlying data X_1 and X_2 follow independent exponential distributions with different parameters.
}
