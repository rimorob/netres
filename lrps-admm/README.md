# lrps-admm: An (accelerated) Alternating Direction Method of Multipliers (ADMM) implementation of the low-rank plus sparse estimator.

This code has now been made into a standalone R package: https://github.com/benjaminfrot/lrpsadmm .

See here for a tutorial: https://nbviewer.jupyter.org/github/benjaminfrot/lrps-admm/blob/master/LRpSADMMExample.ipynb

We implement the estimator described in [Chandrasekaran et al (2012)](https://projecteuclid.org/euclid.aos/1351602527) .

This R implementation uses an accelerated version with restarts of the ADMM algortihm, as described in Goldstein et al  (ftp://ftp.math.ucla.edu/pub/camreport/cam12-35.pdf) .

For questions, remarks or bug reports: benjamin dot frot at gmail dot com . You can also open an issue.

Fitting such methods is challenging from a numerical point-of-view. If you have examples of datasets for which everything went wrong, I am very interested in seeing them; this could be an opportunity to improve the code.

Finally, note that [this MATLAB implementation](https://github.com/TianpeiLuke/LogdetPPA) implements a different kind of algorithm and provides a good implementation of this estimator.
