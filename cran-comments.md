## R CMD check results

0 errors | 0 warnings | 0 notes

## This submission

This update extends the package’s capabilities for spatial empirical dynamic modeling, with a 
focus on methodological flexibility, computational efficiency, and statistical rigor.

* **New features**:
  * Added variable interaction mechanisms in the spatial logistic map for richer coupled 
    dynamics.
  * Introduced configurable distance metrics in cross mapping to adapt distance measures 
    to data characteristics.
  * Enabled alternative spatial cross-sectional embeddings for more flexible state-space 
    reconstruction.

* **Enhancements**:
  * Improved handling of NaNs in cross mapping to ensure robust inference.
  * Safeguarded transient removal in the spatial logistic map to avoid indexing errors.
  * Adjusted default embedding dimension range in `simplex()` and `smap` to `2:10` for 
    stronger attractor reconstruction.
  * Added multithreading in distance-based computations for faster large-scale analysis.
  * Cross mapping plots now display p-value annotations for clearer interpretation.

* **Breaking changes**:
  * Population density dataset coordinates renamed to `lon` and `lat` for consistency.
  * Switched from mean-of-rho to Fisher z-transform for more rigorous significance and 
    confidence interval estimation.

Overall, this release strengthens methodological scope and robustness in spatial EDM workflows.
