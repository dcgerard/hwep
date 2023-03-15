# hwep 2.0.1

- Updates package for compatibility with newer version of rstan.

# hwep 2.0.0

- Added `rmbayes()` and associated functions for Bayes tests for random mating.
- Multiple restarts in EM algorithm for `rmlike()` to fix issues with local convergence.
- Added `rmbayesgl()` for Bayes test of random mating when using genotype likelihoods to account for genotype uncertainty. This uses the `{rstan}` package as a backend.
- Added `glsim()` for genotype likelihood simulation.
- `qqpvalue()` will now allow the user to return a ggplot object.
- Added `menbayesgl()` for testing Mendelian segregation frequencies from S1 and F1 populations.

# hwep 0.0.2

- New function, `qqpvalue()`: A QQ-plot function for $p$-values on the $-\log_{10}$-scale. 
- New function, `ts_bands()`: Calculates simultaneous confidence bands for uniform quantiles.
- New function, `f1dr()`: Estimates the double reduction rate from an F1 population via maximum likelihood.
- Edit aggregation strategy so that it works a little better to mitigate the effects of the asymptotic approximation.

# hwep 0.0.1

- Initial release.
