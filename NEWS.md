# spEDM 1.3

* Implement the `smap` funtion to enable the selection of the optimal theta parameter (#128).

* Now bidirectional mapping in the `gccm` result uses a `full join` structure when organized on the R side (#118).

* Support for calculating unidirectional mappings in the `gccm` function (#117).

* Relax `gccm` C++ source code `libsizes` minimum value constraint of `E+2` (#109).

* Include an option in `gccm` to determine whether to include the current state when generating embedding vectors (#103).

* Provide a complete `GCCM` workflow for spatial lattice and grid data in the `gccm` vignette (#100).

* Add `simplex` function to support selecting the optimal embedding dimension for variables (#98).

* Provide an R-level API for generating embeddings (#97).

* Support testing causal links in GCCM with different `E` and `k` for cause and effect variables (#96).

* Add thread settings for `gccm` (#94).

* Resolve r crash caused by invalid `E` [#90](https://github.com/stscl/spEDM/pull/90) and `k` [#89](https://github.com/stscl/spEDM/pull/89) parameter settings in `gccm`.

* Fix incorrect Pearson correlation calculation in `C++` code when input contains NA (#83).

* Add `S-maps` cross-prediction support to `gccm` (#81).

# spEDM 1.2

* Encapsulate the `gccm` function using the S4 class (#72).

* Add options for `tau`, `k`, and `progressbar` in `gccm` (#69).

* Add `print` and `plot` s3 methods for `gccm` result (#64).

* Require `sdsfun` package version `0.7.0` or higher (#61).

# spEDM 1.1

* Resolve CRAN auto check issues, no significant API changes.

# spEDM 1.0

* Implementing the `GCCM` method for spatial lattice and grid data using pure `C++11`.
