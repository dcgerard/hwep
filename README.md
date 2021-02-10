
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

-   `hwefit()`: Fit either `hwetetra()`,`rmlike()`, or `hwemom()` across
    many loci. Parallelization is supported through the
    [future](https://cran.r-project.org/package=future) package.
-   `hwetetra()`: Likelihood inference for tetraploids. Includes tests
    for random mating, only small deviations from HWE, and the presence
    of double reduction.
-   `rmlike()`: Likelihood inference for random mating in polyploids.
    This function tests for random mating and estimates gametic
    frequencies given random mating. This function does not assume a
    model for meiosis. You should use `hwetetra()` for tetraploids.
-   `hwemom()`: Generalized method of moments inference for double
    reduction. This function tests for equilibrium given double
    reduction rates and estimates these rates given equilibrium. You
    should use `hwetetra()` for tetraploids.

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

## Citation

``` r
citation("hwep")

To cite hwep in publications use:

Gerard D (2021). "Double reduction estimation and equilibrium tests in
natural autopolyploid populations." _Unpublished Manuscript_.

A BibTeX entry for LaTeX users is

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
