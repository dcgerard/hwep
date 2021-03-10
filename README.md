
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

Inference concerning equilibrium and random mating in autopolyploids.
Methods are available to test for equilibrium and random mating at any
even ploidy level (&gt;2) in the presence of double reduction. For
autopolyploid populations in equilibrium, methods are available to
estimate the degree of double reduction. We also provide functions to
calculate genotype frequencies at equilibrium, or after one or several
rounds of random mating, given rates of double reduction.

The main functions for inference are:

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

-   `hwenodr()`: Implements a likelihood ratio test that tests for
    Hardy-Weinberg equilibrium in autopolyploids given no double
    reduction.

Functions are provided for calculating genotype frequencies for
individuals and gametes:

-   `gsegmat()`: Produces the segregation probabilities for gamete
    dosages given parental dosages and the double reduction rate.

-   `gsegmat_symb()`: Provides a symbolic representation of the output
    of `gsegmat()`.

-   `zsegarray()`: Obtains offspring genotype probabilities given
    parental probabilities, the ploidy of the species, and the
    overdispersion parameter, for all possible parental genotypes.

-   `freqnext()`: Updates the genotype frequencies after one generation
    of random mating.

-   `hwefreq()`: Calculate genotype frequencies at equilibrium.

The bounds on the double reduction rate under the complete equational
segregation model are provided by `drbounds()`.

## Installation

<!-- You can install the released version of hwep from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ```{r, eval = FALSE} -->
<!-- install.packages("hwep") -->
<!-- ``` -->

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
