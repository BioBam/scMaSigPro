## scMaSigPro 0.0.5.9000 (development version)

* Breaking changes
    * `as_scmp()` no longer converts Monocle3 `cell_data_set` objects. The `from = "cds"` route and the interactive Shiny app it opened to pick the root and the branching paths are removed, together with the exported `m3_select_path()` and the internal `shinySelect()` and `align_pseudotime()`. `as_scmp()` now takes any `SingleCellExperiment` (subclasses other than `cell_data_set` used to return `NULL`), and `from` defaults to `"sce"`, its only accepted value.
    * The `interactive`, `align_pseudotime` and `anno_col` arguments of `as_scmp()` are removed. They only applied to the Monocle3 route, where `align_pseudotime = TRUE` never worked. Scripts written for earlier versions, including the README, the vignettes and the code accompanying the 2024 paper, pass `align_pseudotime`; delete that argument from the call.
    * The interactive branch selector will be maintained in a separate package that depends on Monocle3. Until then, the removed code is preserved on the `archive/shiny-monocle3` branch (commit 65b5458).
    * To analyse a Monocle3 trajectory now, extract the pseudotime and a path label for every cell and pass the counts and that cell table to `create_scmp()`. The "Common Root Cells" vignette shows the extraction steps.
    * `showCoeff()`, `showGroupCoeff()`, `showInflu()`, `showParams()`, `showSigProf()`, `showSol()` and `showTS()` no longer have `view` and `return` arguments: they always return their result and never open the data viewer. Delete `view =` and `return =` from calls; to browse a result, use `View(showCoeff(obj))`.
    * `.ScMaSigPro_show()`, the helper behind the `show()` method, is no longer exported. Printing an object is unchanged.
    * `sc.p.vector()` and `sc.t.fit()` run their per-gene fits with BiocParallel. The `parallel` and `n_cores` arguments are replaced by `BPPARAM`, which defaults to `BiocParallel::SerialParam()` and fits one gene after another, as `parallel = FALSE` did. Replace `parallel = TRUE, n_cores = k` with `BPPARAM = BiocParallel::MulticoreParam(k)` on Linux and macOS, or with `BPPARAM = BiocParallel::SnowParam(k)` on any system. Windows, where `parallel = TRUE` fell back to one core with a warning, can now fit in parallel with `SnowParam()`. The results are identical with every backend.
    * `sc.p.vector()` and `sc.t.fit()` no longer show a progress bar; with `verbose = TRUE` each reports the number of genes and the backend in one message. For a per-gene progress bar, pass a param created with `progressbar = TRUE`, such as `BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)`.
    * `plotTrendCluster()` loses its `parallel` argument, which had no effect; delete it from calls.
    * `ScMaSigPro` objects have a new `version` slot, the version of scMaSigPro that created them, and an `updateObject()` method brings objects saved by earlier versions up to date. An object saved with `saveRDS()` or `save()` by scMaSigPro 0.0.5 or earlier has no `version` slot: printing it says so, the steps, `pb_counts()`, the plots, `queryCoeff()` and the `show*()` helpers stop with "update it with updateObject() first", and the setters stop with "slots in class definition but not in object: "version"", until it is updated, for example with `ob <- updateObject(readRDS("ob.rds"))`. Its data and results are kept. Objects binned by 0.0.4 cannot be updated because their `Dense` slot is invalid; `updateObject()` stops with the command that bins their cells again with `sc.squeeze()`, which replaces the 0.0.5 advice to run `sc.squeeze()` on them directly. The bundled `scmp.ob` dataset is regenerated with the slot; its other slots are identical.

* Bug Fixes
    * `sc.t.fit()` with `selection_method = "forward"`, `"two.ways.backward"` or `"two.ways.forward"` failed on R >= 4.3 unless `parallel = TRUE` was set, with the error "'length = <n>' in coercion to 'logical(1)'" raised by its progress-bar code. The three methods now run with the default arguments and, like `"backward"`, reproduce `maSigPro::T.fit()` exactly on the package's Gaussian and negative binomial test data.
    * `sc.t.fit()` with `selection_method = "two.ways.backward"` or `"two.ways.forward"` could run forever on some genes, because removing and re-adding terms can cycle; maSigPro's `two.ways.stepback()` and `two.ways.stepfor()` do the same. The selection now stops with the current model when it comes back to an earlier state. Where the copied maSigPro code stopped with an error, the selection now returns a model: tied p-values when a term is removed, a design with at most one estimable term, column names that are not syntactic (forward and two-ways methods), and a two-ways forward step that removes every term (intercept-only model). Every other gene gets the same model as before.
    * The validity methods of `ParameterConfig`, `MatrixDesign`, `VariableProfiles`, `Estimates` and `Significant`, the classes of the slots of an `ScMaSigPro` object, return a description of every problem instead of stopping at the first one, as S4 validity methods should; those of `Estimates` and `Significant` checked nothing. They now check what the workflow relies on: one value in each `Parameters` slot and a `family` object in `distribution`, design matrices and offsets with one row per bin, adjusted p-values named like the p-values, t-score and path-coefficient matrices with the genes of the coefficient matrix, and gene lists made of gene names. Every object the workflow creates passes these checks, so results are unchanged.
    * `validObject()` on an `ScMaSigPro` object reports the problems of every slot, each prefixed with the slot name, instead of stopping with the error of the first invalid slot. `validObject(new("ScMaSigPro"))` now returns `TRUE`; the empty `Sparse` and `Dense` slots of a new object used to be invalid.
    * `sc.set.poly()` failed with "incorrect number of dimensions" when there were several paths and only one bin had the earliest binned pseudotime, as with `create_scmp(use_as_bin = TRUE)` and a single sample at the first time point. It now builds the same design as for any other input.
    * The setters of the `Sparse` and `Dense` slots check their value instead of storing it as given, which could leave an invalid object or metadata that no longer matched the counts. `cSparse<-` (and `cSparse()` with `value`) needs one row per cell, named after the cells and in their order; a table without rownames is taken in the order of the cells and no longer drops the cell names. `eSparse<-` and `eDense<-` need a matrix with the dimensions of the slot whose row and column names, when present, are those of the slot in the same order; missing names are filled in, so `pb_counts(assay_name =)` on such an assay no longer returns zero counts. `cDense<-` still accepts a table with another number of bins, as when bins are rebuilt by hand in the "Common Root Cells" vignette, but it now removes the counts that no longer fit, with a message, so the object stays valid until `pb_counts()` computes the counts of the new bins.
    * The setters (`cSparse<-`, `cDense<-`, `eSparse<-`, `eDense<-`, `pathAssign<-`, `predictors<-`, and `cSparse()`, `cDense()`, `pathAssign()` and `predictors()` with `value`) check the object they return with `validObject()`. `pathAssign<-` and `predictors<-` used to store any matrix, which left an invalid object when its rows were not the bins of the other design matrix (one per bin, with the same names) or when the predictor matrix had another number of columns than the design has terms; such a matrix is now an error and the object is unchanged. Replacing a design matrix by one with the same rows and columns works as before.
    * `as_scmp()`, `create_scmp()`, `sc.set.poly()`, `sc.p.vector()`, `sc.t.fit()`, `sc.filter()` and `sc.cluster.trend()` check the object they return with `validObject()`, as `sc.squeeze()`, `sc.restruct()` and `pb_counts()` already did. Arguments that used to be stored as invalid parameters are now errors naming the parameter: two values for `p_value` in `sc.p.vector()`, `log_offset = NA` with `offset = FALSE` in `sc.p.vector()` and `sc.t.fit()`, and a missing or repeated `root_label` or `path_prefix` in `as_scmp()`. An object broken by hand stops at the next step instead of passing through. Results on valid input are unchanged.
    * The steps, `pb_counts()`, the plots, `queryCoeff()` and the `show*()` helpers check the object they are given before doing any work. An object saved by an earlier version of scMaSigPro stops at once with "update it with updateObject() first" instead of at the end of a step. A `Dense` slot without pseudo-bulk counts for its bins, as after changing the bins with `cDense<-`, asks for `pb_counts()` instead of failing in `assay()` or with "subscript out of bounds". `sc.p.vector()` before `sc.set.poly()`, `sc.t.fit()` when `sc.p.vector()` selected no gene, and `sc.cluster.trend()` before `sc.filter()` name the step to run instead of failing with "'min_na' is too high", "'x' must be an array of at least two dimensions" or an error from `hclust()`. `sc.cluster.trend()`, `plotTrend()`, `plotTrendCluster()` and `plotDiagnostics()` now reject an object of another class with the message of the other functions, which `sc.squeeze()`, `sc.restruct()`, `pb_counts()`, `plotBinTile()` and `plotIntersect()` now share.
    * Arguments that take one of a few values must be one of them, given once, and the column names and labels given to `as_scmp()`, `create_scmp()`, `sc.squeeze()`, `sc.set.poly()`, `sc.restruct()` and `pb_counts()` must be single strings. Two or more values used to pass the check and then fail with "the condition has length > 1" or "Unable to estimate bin size", or to be rejected by `validObject()` after the work was done; they now stop at once with the argument's message. The choices are checked before any work: `bin_method` and `aggregate` of `sc.squeeze()`, `aggregate` of `sc.restruct()` and `pb_counts()`, `selection_method` of `sc.t.fit()` (whose error no longer names the function), `vars` and, with `vars = "groups"`, `intercept` of `sc.filter()` (an unknown `intercept` used to fail with "object 'select_cols' not found"), `geneSet`, `cluster_by`, `cluster_method`, `fill_na` and `use_dim` of `sc.cluster.trend()`, `summary_mode` and, with `logs = TRUE`, `logType` of `plotTrend()`, `model` of `plotDiagnostics()`, and `query` and `change` of `queryCoeff()`.
    * `sc.filter()` with its default `vars` returns the `"all"` gene list instead of failing with "the condition has length > 1". `sc.squeeze()` with `additional_params = NULL` or `list()`, and `as_scmp()` with `additional_params = NULL` or without `labels_exist`, use the defaults instead of failing with "argument is of length zero".
    * `sc.cluster.trend()` checks that the gene set is large enough: at least two genes, as many as `k` for `"hclust"` and more than `k` for `"kmeans"` and `"Mclust"`. A smaller set used to fail inside `hclust()`, `cutree()`, `kmeans()` or `Mclust()`, and a set of one gene clustered with `"kmeans"` or `"Mclust"` returned clusters of the coefficients or bins instead of genes.

* Dependency Changes
    * `shiny`, `plotly`, `igraph` and `scales` are no longer dependencies, so installing scMaSigPro pulls in 55 non-base packages instead of 98.
    * `assertthat` is no longer a dependency. Argument checks use an internal helper and stop with the same messages. A check given `NA` or the wrong number of values (for example two values for an argument that takes one) now reports its own message instead of assertthat's generic "assert_that: ..." text.
    * `stringr` is no longer a dependency; string handling uses base R and gives the same results. Without `assertthat` and `stringr`, installing scMaSigPro pulls in 52 non-base packages instead of 55. `sc.restruct()` still treats `link_sep` as a regular expression, so a separator containing `|` must be escaped (`link_sep = "\\|"`).
    * `e1071`, `entropy` and `MatrixGenerics` are no longer dependencies. Binning computes the skewness (Doane rule) and the equal-width bin counts in base R with the same results, and row sums and means use base R or `Matrix`. `Matrix`, which was always installed with `SingleCellExperiment`, is now declared in Imports.
    * `dplyr` is no longer a dependency, and `rlang` moves from Imports to Suggests because only the tests use it. `sc.squeeze()`, `plotTrend()` and `plotTrendCluster()` group and summarise with base R and give the same bins, counts and plot data; the data of the `lines` layer of `plotTrend()` is a data frame instead of a tibble.
    * `magrittr` is no longer a dependency: the package code no longer uses the `%>%` pipe, which scMaSigPro never exported.
    * `mclust` moves from Imports to Suggests. `sc.cluster.trend(cluster_method = "Mclust")` and the plot of `plotIntersect()` check that `mclust` and `UpSetR` are installed and say how to install them; `plotIntersect(return = TRUE)` does not need `UpSetR`.
    * `ggpubr` is no longer suggested: the vignettes arrange plots with `patchwork`.
    * The ggplot2 functions in use are imported one by one instead of the whole package, and the `dgCMatrix` class is imported from `Matrix`. With these changes installing scMaSigPro pulls in 41 non-base packages instead of 52, and Imports lists 12 packages instead of 18.
    * `parallel` and `parallelly` are replaced by `BiocParallel`, and `utils` is no longer imported. Installing scMaSigPro pulls in 48 non-base packages instead of 41.
    * `maSigPro` moves from Imports to Suggests: only the tests that compare scMaSigPro with maSigPro and the "Bridging maSigPro Analysis" vignette use it, and the vignette shows its code without running it when maSigPro is not installed. Installing scMaSigPro pulls in 44 non-base packages instead of 48 (`maSigPro`, `mclust`, `venn` and `admisc` are no longer installed with it).

* Minor changes
    * Error and progress messages are corrected: misspellings ("doesn't exit", "Overwritting", "onservations", "selction", "select on of"), the class name (`ScMaSigPro`, not `scMaSigPro`), and pointers to functions that do not exist (`plotIntersect()` now points to `sc.filter()`, `showPoly()` to `sc.set.poly()`). `plotDiagnostics()` reports an invalid `model` instead of a missing gene, and messages no longer contain line breaks and indentation copied from the source code. Code that matches the old wording, such as the "Overwritting assay" warning of `eSparse<-` and `eDense<-`, must be updated.
    * `sc.set.poly()` stops with "'poly_degree' must be a number of at least 1" when `poly_degree` is below 1, missing or not a single number. Such values gave a degree-1 design recorded as degree 0 with several paths, and an error in `sc.t.fit()` with one path. A fractional degree still uses its whole part, now also with one path, where a value between 1 and 2 gave a quadratic design recorded as degree 1.
    * `as_scmp()` passes its input through `updateObject()`, so a `SingleCellExperiment` saved by an older version of SingleCellExperiment, such as the bundled `splat.sim` and `multi.lin.sce`, is stored up to date with `labels_exist = FALSE` as it already was with `labels_exist = TRUE`, and `updateObject()` on the new object has nothing to do. Counts and cell metadata are unchanged.

* Internal
    * Unused code is removed: debugging leftovers, commented-out code, the never-enabled `useWeights` branches of `sc.p.vector()` and `sc.t.fit()`, the unexported and broken `ScMaSigPro()` constructor, and 19 unused `importFrom` directives. Results are unchanged.
    * The size factors that `sc.p.vector()` and `sc.t.fit()` use as offsets are computed by a new internal function written from the published median-of-ratios method (Anders and Huber 2010, Genome Biology 11:R106) instead of by a copy of DESeq2's `estimateSizeFactorsForMatrix()`. The copy was LGPL (>= 3) code without credit in a GPL (>= 2) package. No DESeq2 code remains, so the licence stays GPL (>= 2) and no DESeq2 author needs to be added to Authors@R: the function implements the method as published and copies no code. Offsets and results are identical. When every gene has a missing count, `sc.p.vector()` now stops with "Cannot compute size factors" instead of an unrelated error from the model fit.
    * `sc.t.fit()` uses one per-gene function for its four selection methods instead of four copies that differed only in the selection function, and the internal `get_os()`, used only to keep Windows sequential, is removed.
    * The four stepwise selection methods of `sc.t.fit()` are implemented by new internal functions written from the published maSigPro method (Conesa et al. 2006, Bioinformatics 22:1096) instead of by a modified copy of maSigPro's `stepback()`, `stepfor()`, `two.ways.stepback()` and `two.ways.stepfor()`, and the internal `extract_fitting()` no longer calls `maSigPro::position()`. Terms are handled by column position, so the internal `clean_string()` helper is removed. On every design where the copied functions returned a model, the new ones fit the same sequence of models and return an identical fit; the maSigPro oracle tests and the golden results are unchanged.
    * `sc.set.poly()` builds the polynomial design, and `sc.t.fit()` and `plotTrend()` the coefficients of each path, with new internal functions written from the published maSigPro model (Conesa et al. 2006, Bioinformatics 22:1096) instead of maSigPro's `make.design.matrix()` and `reg.coeffs()`. For every whole degree of at least 1 where maSigPro returned a design, the new function returns an identical one, and the path coefficients are identical; the maSigPro oracle tests, the golden results and the plots are unchanged. scMaSigPro no longer calls any maSigPro function.
    * The package reads and writes the `Sparse` and `Dense` slots with the `SingleCellExperiment` accessors `colData()`, `assay()`, `assayNames()` and their replacement forms instead of reaching into the internal slots of those objects, so it no longer depends on how SummarizedExperiment stores its data. Results are unchanged.

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
