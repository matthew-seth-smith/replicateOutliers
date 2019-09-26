[![Build Status](https://travis-ci.org/matthew-seth-smith/replicateOutliers.svg?branch=master)](https://travis-ci.org/matthew-seth-smith/replicateOutliers)

replicateOutliers
===

This package provides functions for outlier detection in the setting of
replicated data. Assuming we have multiple independent measurements (called
"replicates") of each data point, we can find which sets of measurements are so
different from each other that we can call them outliers. We use the absolute
difference
<img src="https://latex.codecogs.com/gif.latex?\Delta"/>
and coefficient of variation
<img src="https://latex.codecogs.com/gif.latex?Z"/>
(Zeta, a kind of relative difference)
for each pair of replicates to make out determinations.

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
methods.

If Travis-CI has an error before completing the initial build, it is because R
for Travis does not suppor the `tcltk` package and I cannot get it to work.
