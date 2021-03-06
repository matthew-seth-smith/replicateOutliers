% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helper_Functions.R
\name{f_gg}
\alias{f_gg}
\title{Joint PDF Under GG Assumption}
\usage{
f_gg(d, z, a_1, b_1, c_1, a_2, b_2, c_2)
}
\arguments{
\item{d}{The absolute difference between two replicates}

\item{z}{The coefficient of variation (CV) between two replicates}

\item{a_1}{The parameter alpha > 0 in the GG distribution for X_1}

\item{b_1}{The parameter beta > 0 in the GG distribution for X_1}

\item{c_1}{The parameter c > 0 in the GG distribution for X_1}

\item{a_2}{The parameter alpha > 0 in the GG distribution for X_2}

\item{b_2}{The parameter beta > 0 in the GG distribution for X_2}

\item{c_2}{The parameter c > 0 in the GG distribution for X_2}
}
\value{
The joint PDF of observing Delta = d and Zeta = z, to be integrated to get useful probabilities
}
\description{
This function is the joint CDF of Delta and Zeta assuming the underlying data X_1 and X_2 follow independent (and possibly identically distributed) GG distributions. We use the parameterization used by the R package flexsurv, given by Kotz and Johnson (1970).
}
