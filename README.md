
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hwep

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN
status](https://www.r-pkg.org/badges/version/hwep)](https://CRAN.R-project.org/package=hwep)
[![R-CMD-check](https://github.com/dcgerard/hwep/workflows/R-CMD-check/badge.svg)](https://github.com/dcgerard/hwep/actions)
[![codecov](https://codecov.io/gh/dcgerard/hwep/branch/main/graph/badge.svg?token=X6QJRSQBXQ)](https://codecov.io/gh/dcgerard/hwep)
<!-- badges: end -->

Inference concerning Hardy-Weinberg equilibrium (HWE) in polyploids.
Methods are available to test for HWE at any ploidy level (&gt;2) in the
presence of double reduction. For autopolyploid populations in HWE,
methods are available to estimate the degree of double reduction. We
also provide functions to calculate genotype frequencies at equilibrium
given rates of double reduction.

The main functions are:

-   `hwefit()`: Fit either `hwelike()`,`rmlike()`, `hweustat()`, or
    `hwenodr()` across many loci. Parallelization is supported through
    the [future](https://cran.r-project.org/package=future) package.
-   `hwelike()`: Likelihood inference for equilibrium. This function
    estimates the rate of double reduction given equilibrium, and tests
    for at most small deviations from equilibrium.
-   `rmlike()`: Likelihood inference for random mating in polyploids.
    This function tests for random mating and estimates gametic
    frequencies given random mating. This function does not assume a
    model for meiosis.
-   `hweustat()`: U-statistic approach for equilibrium and double
    reduction. This function tests for equilibrium given double
    reduction rates and estimates these rates given equilibrium.
-   `hwenodr()`: Implements a likelihood ratio test that tests for HWE
    in autopolyploids given no double reduction.

## Installation

You can install the released version of hwep from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hwep")
```

You can install the development version from
[GitHub](https://github.com/dcgerard/hwep) with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/hwep")
```

## Citation

To cite hwep in publications use:

> Gerard D (2021). “Double reduction estimation and equilibrium tests in
> natural autopolyploid populations.” *Unpublished Manuscript*.

A BibTeX entry for LaTeX users is

``` tex
@Article{,
  title = {Double reduction estimation and equilibrium tests in natural autopolyploid populations},
  author = {David Gerard},
  journal = {Unpublished Manuscript},
  year = {2021},
}
```

## Code of Conduct

Please note that the hwep project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
