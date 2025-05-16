## R CMD check results

0 errors | 0 warnings | 1 note

* These words are correctly spelled:
  - "Gao": a surname of a cited researcher.
  - "granger": refers to "Granger causality", a well-established concept in econometrics.

## This submission

This release includes several new features, improvements, and bug fixes:

* **New R APIs** for false nearest neighbours, spatial granger causality test, and geographical cross mapping cardinality.

* **New vignettes** covering spatial causality methods, SSR for spatial cross-sectional data, and an overview of the package structure.

* **Enhancements** to vignette performance, documentation clarity, and plotting flexibility:
  
  * The S3 plotting method for cross mapping results now supports custom legend text and colors.

* **Breaking changes**:

  * `lib`/`pred` options now exclude `NA` by default.
  
  * Internal example datasets have been refined for consistency with updated method interfaces.
  
* **Bug fixes** addressing C++/R index mismatches and matrix input issues.

This is a regular feature update aimed at expanding modeling support for spatial dynamic empirical modeling workflows and improving overall usability and documentation.
