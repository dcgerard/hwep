
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hwep <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![NSF-2132247](https://img.shields.io/badge/NSF-2132247-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=2132247)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](http://cranlogs.r-pkg.org/badges/grand-total/hwep)](https://cran.r-project.org/package=hwep)
[![CRAN
status](https://www.r-pkg.org/badges/version/hwep)](https://CRAN.R-project.org/package=hwep)
[![R-CMD-check](https://github.com/dcgerard/hwep/workflows/R-CMD-check/badge.svg)](https://github.com/dcgerard/hwep/actions)
[![codecov](https://codecov.io/gh/dcgerard/hwep/branch/main/graph/badge.svg?token=X6QJRSQBXQ)](https://app.codecov.io/gh/dcgerard/hwep)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

Inference concerning equilibrium and random mating in autopolyploids.
Methods are available to test for equilibrium and random mating at any
even ploidy level (\>2) in the presence of double reduction at biallelic
loci. For autopolyploid populations in equilibrium, methods are
available to estimate the degree of double reduction. We also provide
functions to calculate genotype frequencies at equilibrium, or after one
or several rounds of random mating, given rates of double reduction. For
details of these methods, see Gerard (2021).

The main functions for inference are:

-   `hwefit()`: Fit either `hwelike()`,`rmlike()`, `hweustat()`,
    `hwenodr()`, or `hweboot()` across many loci. Parallelization is
    supported through the
    [future](https://cran.r-project.org/package=future) package.

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

-   `hweboot()`: Implements a bootstrap approach to test for equilibrium
    which is more appropriate for small samples and uncertain genotypes.

-   `rmbayes()`: Implements a Bayesian test for random mating in
    autopolyploids for any ploidy level.

-   `rmbayesgl()`: Bayesian test for random mating, accounting for
    genotype uncertainty using genotype likelihoods.

-   `menbayesgl()`: Bayesian test for Mendelian segregation frequencies
    in S1 or F1 populations using genotype likelihoods.

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

Functions for evaluating the uniformity of p-values are provided in
`ts_bands()` and `qqpvalue()`.

## Installation

You can install the released version of hwep from
[CRAN](https://cran.r-project.org/package=hwep) with:

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

> Gerard D (2022). “Double reduction estimation and equilibrium tests in
> natural autopolyploid populations.” *Biometrics* In press.
> [doi:10.1111/biom.13722](https://doi.org/10.1111/biom.13722).

A BibTeX entry for LaTeX users is

``` tex
@Article{,
  title = {Double reduction estimation and equilibrium tests in natural autopolyploid populations},
  author = {David Gerard},
  journal = {Biometrics},
  year = {2022},
  doi = {10.1111/biom.13722},
  volume = {In press},
}
```

If you use `rmbayes()`, `rmbayesgl()`, or `menbayeslg()`, then please
also cite

> Gerard D (2022). “Bayesian tests for random mating in autopolyploids.”
> *bioRxiv*.
> [doi:10.1101/2022.08.11.503635](https://doi.org/10.1101/2022.08.11.503635).

A BibTeX entry for LaTeX users is

``` tex
@article{,
    author = {Gerard, David},
    title = {Bayesian Tests for Random Mating in Autopolyploids},
    year = {2022},
    doi = {10.1101/2022.08.11.503635},
    publisher = {Cold Spring Harbor Laboratory},
    journal = {bioRxiv}
}
```

## Acknowledgments

This material is based upon work supported by the National Science
Foundation under Grant
No. [2132247](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2132247).
The opinions, findings, and conclusions or recommendations expressed are
those of the author and do not necessarily reflect the views of the
National Science Foundation.

## Code of Conduct

Please note that the hwep project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
