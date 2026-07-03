# spEDM

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![CRAN Release](https://www.r-pkg.org/badges/last-release/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![CRAN Checks](https://badges.cranchecks.info/worst/spEDM.svg)](https://cran.r-project.org/web/checks/check_results_spEDM.html)
[![Downloads_all](https://badgen.net/cran/dt/spEDM?color=orange)](https://CRAN.R-project.org/package=spEDM)
[![Downloads_month](https://cranlogs.r-pkg.org/badges/spEDM)](https://CRAN.R-project.org/package=spEDM)
[![License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-20b2aa.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/stscl/spEDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stscl/spEDM/actions/workflows/R-CMD-check.yaml)
[![R-universe](https://stscl.r-universe.dev/badges/spEDM?color=cyan)](https://stscl.r-universe.dev/spEDM)
[![IJGIS](https://img.shields.io/badge/IJGIS-10.1080%2F13658816.2026.2687121-f25b3e?logo=doi&style=flat)](https://doi.org/10.1080/13658816.2026.2687121)

<!-- badges: end -->

<a href="https://stscl.github.io/spEDM/"><img src="man/figures/spEDM.png" align="right" hspace="5" vspace="0" width="15%" alt="spEDM website: https://stscl.github.io/spEDM/"/></a>

***Sp**atial **E**mpirical **D**ynamic **M**odeling*

*spEDM* is an R package for spatial causal discovery. It extends Empirical Dynamic Modeling (EDM) from time series to spatial cross-sectional data, provides seamless support for vector and raster spatial data via tight integration with the [*sf*](https://CRAN.R-project.org/package=sf) and [*terra*](https://CRAN.R-project.org/package=terra) packages, and enables data-driven causal inference from spatial snapshots.

> *Refer to the package documentation <https://stscl.github.io/spEDM/> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=spEDM) with:

``` r
install.packages("spEDM", dependencies = TRUE)
```

- Install binary version from [R-universe](https://stscl.r-universe.dev/spEDM) with:

``` r
install.packages("spEDM",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dependencies = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/spEDM) with:

``` r
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}
pak::pak("stscl/spEDM", dependencies = TRUE)
```

## CITATION

Please cite **[spEDM][1]** as:

```
Lyu, W., Dai, S., Song, Y., Zhao, W., Yi, W., Xiao, Y., Jia, N., 2026. Measuring causal strengths from spatial cross-sectional data with geographical cross mapping cardinality. International Journal of Geographical Information Science 1–23. https://doi.org/10.1080/13658816.2026.2687121
```

A BibTeX entry for LaTeX users is:

``` bib
@article{lyu2026gcmc, 
    title = {Measuring causal strengths from spatial cross-sectional data with geographical cross mapping cardinality}, 
    ISSN = {1362-3087}, 
    DOI = {10.1080/13658816.2026.2687121}, 
    journal = {International Journal of Geographical Information Science}, 
    publisher = {Informa UK Limited}, 
    author = {Lyu, Wenbo and Dai, Shaoqing and Song, Yongze and Zhao, Wufan and Yi, Wen and Xiao, Yumiao and Jia, Nan}, 
    year = {2026}, 
    month = {June}, 
    pages = {1–23} 
}
```

## Reference

Lyu, W., Dai, S., Song, Y., Zhao, W., Yi, W., Xiao, Y., Jia, N., 2026. Measuring causal strengths from spatial cross-sectional data with geographical cross mapping cardinality. International Journal of Geographical Information Science 1–23. [https://doi.org/10.1080/13658816.2026.2687121][1].

Lyu, W., Lei, Y., Yi, W., Song, Y., Li, X., Dai, S., Qin, Y., Zhao, W., 2026. Causal discovery in urban data with temporal empirical dynamic modeling: The R package tEDM. Computers, Environment and Urban Systems 127, 102435. [https://doi.org/10.1016/j.compenvurbsys.2026.102435][2].

Gao, B., Yang, J., Chen, Z., Sugihara, G., Li, M., Stein, A., Kwan, M.-P., Wang, J., 2023. Causal inference from cross-sectional earth system data with geographical convergent cross mapping. Nature Communications 14. [https://doi.org/10.1038/s41467-023-41619-6][3].

Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., Munch, S., 2012. Detecting Causality in Complex Ecosystems. Science 338, 496–500. [https://doi.org/10.1126/science.1227079][5].

&nbsp; 

[1]: https://doi.org/10.1080/13658816.2026.2687121
[2]: https://doi.org/10.1016/j.compenvurbsys.2026.102435
[3]: https://doi.org/10.1038/s41467-023-41619-6
[4]: https://doi.org/10.1126/science.1227079
[5]: https://doi.org/10.1126/science.1227079
