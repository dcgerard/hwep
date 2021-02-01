
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hwep

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN
status](https://www.r-pkg.org/badges/version/hwep)](https://CRAN.R-project.org/package=hwep)
[![R-CMD-check](https://github.com/dcgerard/hwep/workflows/R-CMD-check/badge.svg)](https://github.com/dcgerard/hwep/actions)
[![codecov](https://codecov.io/gh/dcgerard/hwep/branch/main/graph/badge.svg?token=X6QJRSQBXQ)](https://codecov.io/gh/dcgerard/hwep)
<!-- badges: end -->

Inference concerning Hardy-Weinberg equilibrium (HWE) in polyploids.
Methods are available to test for HWE at any ploidy level in the
presence of double reduction. For autopolyploid populations in HWE,
methods are available to estimate the degree of double reduction. We
also provide functions to calculate genotype frequencies at equilibrium
given rates of double reduction.

The main functions are:

-   `hwelike()`: Likelihood inference for HWE in polyploids. This
    function tests for HWE and estimates gametic frequencies given HWE.
    This function does not assume a model for meiosis.
-   `hwemom()`: Generalized method of moments inference for double
    reduction. This function tests for equilibrium given double
    reduction rates and estimates these rates given equilibrium. Note
    that for ploidies greater than four, not all gametic frequencies are
    possible given this model for meiosis. That is, the test might
    reject not because the population is not in HWE, but because the
    model for meiosis is incorrect. So `hwelike()` is more general if
    you just want to test for HWE and do not care about double
    reduction.

## Installation

You can install the released version of hwep from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hwep")
```

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/hwep")
```

## Code of Conduct

Please note that the hwep project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
