# spEDM 1.3

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
