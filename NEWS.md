## scMaSigPro 0.0.4.9000 (development version)

* Documentation
    * Vignettes no longer download data or images from the web. The "Common Root Cells" vignette uses the bundled `multi.lin.sce` dataset and includes the Monocle3 extraction steps as an unevaluated recipe; the Quick Start vignette no longer writes into `data/` during the build.
    * Installation instructions use `BiocManager::install("BioBam/scMaSigPro")` instead of pinning Bioconductor 3.14.
    * Scripts that regenerate the bundled datasets live in `inst/scripts/`; the dataset help pages point to them instead of embedding the code.

* Infrastructure
    * `.Rbuildignore` is now tracked, so fresh clones and CI build the same package as the maintainer.
    * GitHub Actions rebuilt on the current r-lib templates: R CMD check (plus a weekly run), test coverage via codecov, lint (advisory), pkgdown, and the R-universe build that mirrors Bioconductor's submission checks including BiocCheck. Workflows run on `main` and `devel`.
    * DESCRIPTION: title in title case, Description cites the paper, ORCIDs for all four authors (the previous `OCRID` typo meant none were recognised), `LazyData` removed, roxygen markdown enabled.

## scMaSigPro 0.0.4 (2024-07-30)

* Dependency Changes
    * Removed hard dependency on `RColorConesa` and `ComplexUpset`.
    * Removed `assertthat`, `stringr` from namespace.
    * Removed `maSigPro::position` from namespace.

* Plotting Updates
    * Updated plot functions to use a custom palette.
    * Updated `plotTrendCluster` function for distinct layers.
    * Enhanced `plotIntersect()` with return capabilities.
    * Enhanced `plotTrend` with curves, lines, and points.
    * Updated ordering by frequency in `UpsetR`.

* Function Updates
    * Added `clean_string()` internal function.
    * Added `sc.restruct()` function.
    * Exported `pb_helpers()`.

* Documentation Updates
    * General code styling.
    * Updated README with citation information.
    * Single comprehensive data documentation.
    * Added a new vignette.

* Bug Fixes
    * Removed `patchwork` call.

## scMaSigPro 0.0.3 (2024-03-3)

* Bug fixes
    * `eSparse()` and `eDense()` indexes
    * Removed the ComplexUpset error due to ggplot 3.5.0.

* Remove Package Dependence
    * SummarizedExperiment
    * ComplexUpset
    
* Package Suggestions Added
    * Patchwork
    * ComplexUpset
    * UpSetR

* Updated Test cases
    * Functionality of `eSparse()` and `eDense()`
    * Testing `sc.squeeze()` against manual pseudo-bulking

## scMaSigPro 0.0.1 (2023-12-20)

* Initial Github release.
