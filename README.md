
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spEDM <img src="man/figures/logo.png" align="right" height="139" alt="https://stscl.github.io/spEDM/">

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![CRAN
Release](https://www.r-pkg.org/badges/last-release/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![CRAN
Checks](https://badges.cranchecks.info/worst/spEDM.svg)](https://cran.r-project.org/web/checks/check_results_spEDM.html)
[![Downloads_all](https://badgen.net/cran/dt/spEDM?color=orange)](https://CRAN.R-project.org/package=spEDM)
[![Downloads_month](https://cranlogs.r-pkg.org/badges/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![R-universe](https://stscl.r-universe.dev/badges/spEDM?color=cyan)](https://stscl.r-universe.dev/spEDM)
[![R-CMD-check](https://github.com/stscl/spEDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stscl/spEDM/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-20b2aa.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

<!-- badges: end -->

***Sp**atial **E**mpirical **D**ynamic **M**odeling*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=spEDM) with:

``` r
install.packages("spEDM", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/spEDM) with:

``` r
install.packages("spEDM",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/spEDM)
  with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/spEDM",
                         build_vignettes = TRUE,
                         dep = TRUE)
```
