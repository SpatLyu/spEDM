# spEDM 1.7

### new

* Provide R-level API and vignette for spatial logistic map ([#609](https://github.com/stscl/spEDM/pull/609),[#612](https://github.com/stscl/spEDM/pull/612)).

### enhancements

* Prevent duplicate generic registrations and simplify method definitions (#627).

* Unify parameter descriptions to lowercase (#570).

* Introduce confidence interval ribbon support in plot method for `gccm` results (#550).

### breaking changes

* Symbolization `C++` functions now compute medians from lib subset only (#599).

* Users must now use `detrend` instead of `trend.rm` (#559).

* Enable `column` parameter in `simplex()` and `smap()` and rename `columns` parameter to `column` in `multiview()` (#565).

### bug fixes

* Fix bug in plot method for `gccm` results where legend labels did not correctly match the corresponding line colors (#552).

# spEDM 1.6

### new

* Enact `fnn` R API support for false nearest neighbours method (#512).

* Integrate R API and vignette for geographical cross mapping cardinality (#455).

* Provide R-level API and vignette for spatial causality test (#403).

### enhancements

* Enable custom legend texts and colors in plot method for `gccm` results (#535).

* Reduce computational load in vignettes (#476).

* Document overall structure and usage of spEDM in a dedicated vignette (#415).

* Create `SSR` vignette for spatial cross-sectional data state-space reconstruction (#412).

* Include references for algorithms in `spEDM` (#367).

### breaking changes

* Use non-NA spatial units for `lib`/`pred` by default (#499).

* Refine internal case data (#348).

### bug fixes

* Patch memory error caused by mismatch between C++ (0-based) and R (1-based) indexing (#480).

* Fix error from non-matrix input in grid-type handling due to R matrix slicing (#474).

# spEDM 1.5

### new

* Enable `parallel.level` parameter to specify parallel granularity in `gccm` R API (#310).

* Implement the `multiview` function for multiview embedding forecasting (MVE) method (#221).

### enhancements

* Integrate `lib` parameter in `gccm` R API for library units selection (#278).

* Set the default `k` to `E+2` in the `gccm` R API (#261).

* Eliminate redundant computations at the source C++ code level (#233).

* Add `trend.rm` option in the R API for `embedded`, `simplex`, and `smap` methods to align with `gccm` behavior (#191).

* Refactor indexing of lag values and embedding vector generation for spatial lattice ([#186](https://github.com/stscl/spEDM/pull/186),[#184](https://github.com/stscl/spEDM/pull/184)) and grid data ([#183](https://github.com/stscl/spEDM/pull/183),[#181](https://github.com/stscl/spEDM/pull/181)).

### breaking changes

* Default plotting method places the legend in the top-left corner of the plot now (#325).

* Refine `simplex` & `smap` output on the R side (#263).

### bug fixes

* Fix bug in R functions `embedded`, `simplex`, `smap` when input data contains only one attribute column (#246).

# spEDM 1.4

### enhancements

* Improve default spatial neighbors list generation for spatial lattice data with support from the `sdsfun` package (#159).

### breaking changes

* Adjust the behavior of the `tau` parameter in the C++ source code and update the R side API (#154).

# spEDM 1.3

### new

* Implement the `smap` function to enable the selection of the optimal theta parameter (#128).

* Add `simplex` function to support selecting the optimal embedding dimension for variables (#98).

* Provide an R-level API for generating embeddings (#97).

### enhancements

* Now bidirectional mapping in the `gccm` result uses a `full join` structure when organized on the R side (#118).

* Support for calculating unidirectional mappings in the `gccm` function (#117).

* Relax `gccm` C++ source code `libsizes` minimum value constraint of `E+2` (#109).

* Provide a complete `GCCM` workflow for spatial lattice and grid data in the `gccm` vignette (#100).

* Support testing causal links in GCCM with different `E` and `k` for cause and effect variables (#96).

* Add thread settings for `gccm` (#94).

* Add `S-maps` cross-prediction support to `gccm` (#81).

### bug fixes

* Resolve r crash caused by invalid `E` [#90](https://github.com/stscl/spEDM/pull/90) and `k` [#89](https://github.com/stscl/spEDM/pull/89) parameter settings in `gccm`.

* Fix incorrect Pearson correlation calculation in `C++` code when input contains NA (#83).

# spEDM 1.2

### enhancements

* Encapsulate the `gccm` function using the S4 class (#72).

* Add options for `tau`, `k`, and `progressbar` in `gccm` (#69).

* Add `print` and `plot` s3 methods for `gccm` result (#64).

### bug fixes

* Fix the bug where the `gccm` function returns empty results when input grid data contains NA values (#61).

# spEDM 1.1

### bug fixes

* Resolve CRAN auto check issues, no significant API changes.

# spEDM 1.0

### new

* Implementing the `GCCM` method for spatial lattice and grid data using modern C++.
