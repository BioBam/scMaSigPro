## scMaSigPro 0.0.5.9000 (development version)

* Breaking changes
    * `as_scmp()` no longer converts Monocle3 `cell_data_set` objects. The `from = "cds"` route and the interactive Shiny app it opened to pick the root and the branching paths are removed, together with the exported `m3_select_path()` and the internal `shinySelect()` and `align_pseudotime()`. `as_scmp()` now takes any `SingleCellExperiment` (subclasses other than `cell_data_set` used to return `NULL`), and `from` defaults to `"sce"`, its only accepted value.
    * The `interactive`, `align_pseudotime` and `anno_col` arguments of `as_scmp()` are removed. They only applied to the Monocle3 route, where `align_pseudotime = TRUE` never worked. Scripts written for earlier versions, including the README, the vignettes and the code accompanying the 2024 paper, pass `align_pseudotime`; delete that argument from the call.
    * The interactive branch selector will be maintained in a separate package that depends on Monocle3. Until then, the removed code is preserved on the `archive/shiny-monocle3` branch (commit 65b5458).
    * To analyse a Monocle3 trajectory now, extract the pseudotime and a path label for every cell and pass the counts and that cell table to `create_scmp()`. The "Common Root Cells" vignette shows the extraction steps.
    * `showCoeff()`, `showGroupCoeff()`, `showInflu()`, `showParams()`, `showSigProf()`, `showSol()` and `showTS()` no longer have `view` and `return` arguments: they always return their result and never open the data viewer. Delete `view =` and `return =` from calls; to browse a result, use `View(showCoeff(obj))`.
    * `.ScMaSigPro_show()`, the helper behind the `show()` method, is no longer exported. Printing an object is unchanged.

* Dependency Changes
    * `shiny`, `plotly`, `igraph` and `scales` are no longer dependencies, so installing scMaSigPro pulls in 55 non-base packages instead of 98.

* Minor changes
    * Error and progress messages are corrected: misspellings ("doesn't exit", "Overwritting", "onservations", "selction", "select on of"), the class name (`ScMaSigPro`, not `scMaSigPro`), and pointers to functions that do not exist (`plotIntersect()` now points to `sc.filter()`, `showPoly()` to `sc.set.poly()`). `plotDiagnostics()` reports an invalid `model` instead of a missing gene, and messages no longer contain line breaks and indentation copied from the source code. Code that matches the old wording, such as the "Overwritting assay" warning of `eSparse<-` and `eDense<-`, must be updated.

* Internal
    * Unused code is removed: debugging leftovers, commented-out code, the never-enabled `useWeights` branches of `sc.p.vector()` and `sc.t.fit()`, the unexported and broken `ScMaSigPro()` constructor, and 19 unused `importFrom` directives. Results are unchanged.

## scMaSigPro 0.0.5

* Bug Fixes
    * `sc.squeeze()` no longer scrambles the cell metadata in the `Sparse` slot. Binning stacked the cells path by path, but the count matrix kept its original column order, so `colData()`, `counts()`, `assay()` and `cSparse()` paired every cell with another cell's counts (on `splat.sim`, `counts(ob@Sparse)[, "Cell2"]` returned the counts of `Cell1`). Pseudo-bulk counts, fits and gene lists were not affected because they match cells by name. Running `sc.squeeze()` again repairs objects binned with 0.0.4.
    * The `Dense` slot is now a valid `SingleCellExperiment`. It used to report 0 genes (`dim(ob@Dense)` was `0 x <bins>`), failed `validObject()`, and `ob@Dense[1:3, ]` errored. `pb_counts()` now builds it from the `bulk.counts` assay and the bin metadata together; `sc.squeeze()`, `sc.restruct()` and `pb_counts()` check the result with `validObject()`.
    * `pb_counts()` names the bins after the rownames of the `Dense` metadata (the `bin_col` column is used only when there are no rownames) and stops on duplicated names. After bins were reordered by hand with `cDense<-` (Common Root Cells vignette), the counts kept the old, duplicated bin labels and `sc.p.vector()` failed with "subscript out of bounds"; such objects now run through the workflow.
    * When `split_bins = TRUE` leaves cells outside every bin (3 of the 1,977 cells of `multi.lin.sce` with the default `bin_method`), `sc.squeeze()` now warns and removes them from the `Sparse` slot. Previously they were silently dropped from the cell metadata only, which left the object invalid. They were never part of any bin, so results do not change.
    * The bundled `scmp.ob` dataset was regenerated: it passes `validObject()` and its results are numerically unchanged.
    * Printing a binned object lists the branching paths in the order they appear in the cell metadata, as before binning.

* Dependency Changes
    * `SummarizedExperiment` (already installed with `SingleCellExperiment`) is declared in Imports.

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
