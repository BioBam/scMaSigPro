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
