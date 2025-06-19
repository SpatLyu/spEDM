## R CMD check results

0 errors | 0 warnings | 0 notes

## This submission

This release introduces new modeling capabilities, streamlined method infrastructure, and refinements to the package interface and visualization support:

* **New R API and vignette** for the spatial logistic map, enabling simulation and analysis of spatial nonlinear dynamics.

* **Enhancements** to method infrastructure, documentation consistency, and plotting functionality:

  * Replaced logical vectors with integer indices for more efficient `lib/pred` indexing in `C++` level forecasting methods.

  * Generic method registration has been made more robust, preventing duplicate registrations and simplifying S4 method management.

  * All parameter descriptions have been standardized to lowercase for improved documentation clarity.

  * Confidence interval ribbons are now supported in the S3 plotting method for cross mapping results.

* **Breaking changes**:

  * Refined randomization strategy in spatial causality test.

  * C++ symbolization functions now compute medians using only the `lib` subset, ensuring consistency with library-based inference.

  * The argument `trend.rm` has been renamed to `detrend` for improved naming consistency.

  * The `column` argument is now uniformly supported in `simplex()`, `smap()`, and `multiview()` (renaming `columns` to `column`).

* **Bug fixes**:

  * Corrected a mismatch between legend labels and line colors in plots of cross mapping results.

This is a regular feature update focused on extending support for spatial nonlinear modeling workflows, improving code robustness, and maintaining interface consistency.
