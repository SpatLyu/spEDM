
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spEDM <img src="man/figures/logo.png" align="right" height="139" alt="https://stscl.github.io/spEDM/">

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-cyan.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN](https://www.r-pkg.org/badges/version/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![DownloadsAll](https://badgen.net/cran/dt/spEDM?color=orange)](https://CRAN.R-project.org/package=spEDM)
[![DownloadsMonth](https://cranlogs.r-pkg.org/badges/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![R-universe](https://stscl.r-universe.dev/badges/spEDM?color=cyan)](https://stscl.r-universe.dev/spEDM)
<!-- badges: end -->

**Spatial Empirical Dynamic Modeling**

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
