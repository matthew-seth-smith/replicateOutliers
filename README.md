[![Build Status](https://travis-ci.org/matthew-seth-smith/replicateOutliers.svg?branch=master)](https://travis-ci.org/matthew-seth-smith/replicateOutliers)

replicateOutliers
===

This package provides functions for outlier detection in the setting of
replicated data. Assuming we have multiple independent measurements (called
"replicates") for each data point, we can find which sets of measurements are so
different from each other that we can call them outliers. We use the absolute
difference
<img src="https://latex.codecogs.com/gif.latex?\Delta"/>
and coefficient of variation
<img src="https://latex.codecogs.com/gif.latex?Z"/>
(Zeta, a kind of relative difference)
for each pair of replicates to make our determinations.

The function `outlier_DZ` returns a numeric identifier for the outlier status
(
<img src="https://latex.codecogs.com/gif.latex?0"/>
for non-outlier,
<img src="https://latex.codecogs.com/gif.latex?1"/>
for large
<img src="https://latex.codecogs.com/gif.latex?\Delta"/>
but not
<img src="https://latex.codecogs.com/gif.latex?Z"/>,
and
<img src="https://latex.codecogs.com/gif.latex?2"/>
for
outlier). The other main functions are `q_exp_joint_DZ`, `q_exp_marg_DZ`,
`q_gg_joint_DZ`, and `q_gg_marg_DZ`. Their outputs are probabilities called
<img src="https://latex.codecogs.com/gif.latex?q"/>-values
that we determine using joint or marginal probability distributions.

Please see {CITE PAPER} for further details and the derivations of these
methods. Also please check the vignette HTML file in the `inst/doc` directory.

Currently, `outlier_DZ` does not work, because it depends on the `extremevalues`
package for the `getOutliersI` function. This package, in turn, depends on the
`tcltk` package, which is now part of base `R`. The package `extremevalues`
should then be able to run without importing `tcltk`, but it does not work. If
you are able to get `extremevalues` to work locally, the code for `outlier_DZ`
is in the vignette explaining this method.

To install this package, please install and load the `devtools` package and use
the command

```{R}
install_github("matthew-seth-smith/replicateOutliers")
```

Travis-CI is throwing errors for this build (as does loading this function on
another computer) because of conflicts between the `gompertz` functions used by
the packages `VGAM` and `flexsurv`. I do not believe this package uses the
Gompertz Distribution, so this should not be a problem.
