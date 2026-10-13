# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`scMaSigPro` is an R package (S4, Bioconductor-style) for differential expression along single-cell trajectories. It pseudo-bulks cells along pseudotime into bins, then runs the maSigPro polynomial-GLM workflow (negative binomial with size-factor offsets) to find genes whose expression profiles differ between branching paths. Published as Srivastava et al., Bioinformatics 2024 (btae443). Version 0.0.4, GPL >= 2, installed from GitHub (not CRAN/Bioconductor).

## Commands

Run from the package root in an R session or via `Rscript -e '...'`. Local toolchain has R 4.6 with devtools, testthat, roxygen2, styler, lintr; `covr` is not installed.

| Task | Command |
|---|---|
| Bioconductor deps | `BiocManager::install(c("SingleCellExperiment","maSigPro","MatrixGenerics","S4Vectors"))` |
| Load for development | `devtools::load_all()` |
| Regenerate `man/` + `NAMESPACE` | `devtools::document()` (or `roxygen2::roxygenise()`) |
| All tests | `devtools::test()` |
| One test file | `devtools::test(filter = "reproduce-masigpro")` or `testthat::test_file("tests/testthat/test-reproduce-masigpro.R")` |
| Full check (what CI runs) | `devtools::check()` or `rcmdcheck::rcmdcheck()` |
| Style | `styler::style_pkg()` |
| Lint | `lintr::lint_package()` |
| Install locally | `devtools::install(build_vignettes = FALSE)` |
| Build vignettes | `devtools::build_vignettes()` (needs BiocStyle, knitr, rmarkdown; two vignettes download objects from metapriyansh.com) |

Notes:
- `NAMESPACE` and `man/` are roxygen output. Never hand-edit them; change roxygen tags in `R/` and run `devtools::document()`.
- `assertthat` and `stringr` are in `Imports` but deliberately not imported into NAMESPACE. Call them as `assertthat::assert_that()` / `stringr::str_*()`. `maSigPro::position()` is likewise called with `::` on purpose.
- `UpSetR` and `ggpubr` are only in `Suggests`, yet `plotIntersect()` calls `UpSetR::upset()` unconditionally.
- `.Rbuildignore` is gitignored, so a fresh clone has none. Expect "non-standard top-level files" NOTEs from `R CMD check` for `renv.lock`, `pkgdown/`, `scMaSigPro.Rproj`.
- `renv.lock` is a snapshot only (R 4.3 / Bioc 3.17). There is no `.Rprofile` activating renv.
- CI (`.github/workflows/`): R-CMD-check matrix on macOS/Windows/Ubuntu release plus Ubuntu devel, for pushes to `main` and PRs to `main`/`develop`; `covr::codecov()` coverage; PR-comment commands `/document` and `/style` (maintainers only) run roxygen/styler and push a commit.
- Tests use testthat edition 3. `test-reproduce-masigpro.R` needs the `maSigPro` package and its bundled `data.abiotic` / `edesign.abiotic` datasets.

## Architecture

### The `ScMaSigPro` S4 object (R/ScMaSigProClass.R)

Every pipeline step takes the object and returns it with one more slot filled. Seven slots:

| Slot | Class | Filled by | Holds |
|---|---|---|---|
| `Sparse` | `SingleCellExperiment` | `as_scmp()` / `create_scmp()` | cell-level `counts` assay + colData (pseudotime, path) |
| `Dense` | `SingleCellExperiment` | `sc.squeeze()` / `sc.restruct()` | bin-level `bulk.counts` assay + bin colData |
| `Design` | `MatrixDesign` | `sc.set.poly()` | `predictor_matrix` (maSigPro `dis`), `assignment_matrix` (`edesign`), `groups.vector`, `offset` |
| `Profile` | `VariableProfiles` | `sc.p.vector()` | per-gene `p_values`, `adj_p_values`, `fdr`, `non_flat` gene names |
| `Estimate` | `Estimates` | `sc.t.fit()` | `significance_matrix` (sol), `coefficient_matrix`, `path_coefficient_matrix`, `t_score_matrix`, `influential` |
| `Significant` | `Significant` | `sc.filter()` / `sc.cluster.trend()` | `genes` (named list) and `clusters` |
| `Parameters` | `ParameterConfig` | every step | colData column names + every fitting parameter used so far |

`@Parameters` is the shared state. Each step records what it used (`poly_degree`, `p_value`, `distribution`, `offset`, `link`, ...) and downstream functions default their arguments from it, e.g. `sc.t.fit(p_value = scmpObj@Parameters@p_value)`. It stores the names of colData columns, not the data: `ptime_col`, `path_col`, `bin_ptime_col` (`scmp_binned_pseudotime`), `bin_col` (`scmp_bin`), `bin_size_col` (`scmp_bin_size`), `bin_mem_col` (`scmp_bin_members`). Look columns up through these slots, never by literal name.

Accessors (R/scmpGenerics.R): `cSparse()` / `cDense()` return colData as a data.frame; `eSparse()` / `eDense()` return an assay as a dense matrix (with the default `i`, they return the first assay regardless of its name); `pathAssign()` / `predictors()` expose the Design matrices. Setters coerce to `dgCMatrix`. The `show*()` helpers in R/show_functions.R read `@Estimate`: `showSol()`, `showCoeff()`, `showTS()`, `showGroupCoeff()`, `showInflu()`, `showSigProf()`, plus `showPoly()` (formula string) and `showParams()`.

### Pipeline

```
as_scmp(sce)   or   create_scmp(counts, cell_data, ptime_col, path_col)
  -> sc.squeeze()        bin pseudotime per path            -> @Dense          R/sc.squeeze.R, R/pb_helpers.R
     -> sc.restruct()    optional: multi-branch, shared root                    R/sc.restruct.R
  -> sc.set.poly()       maSigPro::make.design.matrix        -> @Design         R/sc.set.poly.R
  -> sc.p.vector()       full GLM vs intercept, per gene     -> @Profile        R/sc.p.vector.R
  -> sc.t.fit()          stepwise term selection             -> @Estimate       R/sc.t.fit.R, R/maSigPro_term_selection.R
  -> sc.filter()         R2 / p-value gene lists             -> @Significant    R/sc.filter.R
  -> sc.cluster.trend()  hclust / kmeans / Mclust            -> @Significant    R/sc.cluster.trend.R
  -> plotTrend(), plotTrendCluster(), plotIntersect(), plotBinTile(), plotDiagnostics(), queryCoeff()
```

**Pseudo-bulking (`sc.squeeze`)**: per path, estimates a bin count from the pseudotime vector (`estBinSize()` in R/pb_helpers.R: Sturges, Freedman-Diaconis, Sqrt, Rice, Doane, Scott), discretizes with `entropy::discretize`, then optionally splits oversized bins (`split_bins`), prunes small ones (`prune_bins`), renumbers gaps (`fill_gaps`) and drops trailing bins so paths have equal length (`drop_trails`). Each bin's member cells are stored in Dense colData as a `|`-separated string in `scmp_bin_members`; `pb_counts()` parses that string to sum or mean the Sparse counts into `bulk.counts`. Bin rownames are `<path>_bin_<n>`.

**maSigPro lineage**: `sc.set.poly` dummy-codes paths, derives a `Replicate` column by run-length encoding, and calls `maSigPro::make.design.matrix`. `sc.p.vector` / `sc.t.fit` are re-implementations of `maSigPro::p.vector` / `T.fit` that add a GLM `offset` (DESeq2 median-of-ratios size factors, vendored in R/DESeq2_estimateSizeFactorsForMatrix.R, computed on counts + 1) and default to `MASS::negative.binomial(10)` with a log link. R/maSigPro_term_selection.R holds copies of maSigPro's `stepback` / `stepfor` / `two.ways.*` modified to pass offset and weights. **Invariant**: with `family = gaussian()`, `link = "identity"`, `offset = FALSE`, results must be `identical()` to maSigPro's. `tests/testthat/test-reproduce-masigpro.R` enforces this slot by slot (design matrices, p-values, FDR, sol, coefficients, t-scores, influence, gene lists). Run it after touching any fitting code.

**Trajectory-tool outputs**: the package takes a `SingleCellExperiment` whose `colData` already holds a pseudotime column and a path column (`as_scmp()`), or a count matrix plus a cell table with those two columns (`create_scmp()`). Users extract pseudotime and hard path labels from their trajectory tool first; the Common-Root-Cells vignette shows the Monocle3 recipe. The former interactive Shiny selector and the `from = "cds"` route (`R/shinySelect.R`, `R/m3_select_path.R`, `R/align_pseudotime.R`) were removed in 0.0.5.9000 and are preserved on the `archive/shiny-monocle3` branch (commit 65b5458) for a future companion package. Never reference `monocle3` from `R/`: it is GitHub-only and an undeclared namespace is an R CMD check WARNING.

**Multi-branch with a common root (`sc.restruct`)**: after `sc.squeeze`, duplicates the root path's bins into each leaf path and shifts each leaf's `bin_ptime_col` by the root's bin count, so the design sees full-length paths that share early bins. Added in PR #2 (`multi_lineage`); see `vignettes/Common-Root-Cells.Rmd` and the `multi.lin.sce` dataset.

### Bundled data (`data/`, documented in R/data_doc.R)

- `splat.sim`: 200 cells x 100 genes, splatter bifurcation; pseudotime in `Step`, path in `Group`. Tests and the basic vignette load it via `as_scmp(..., additional_params = list(labels_exist = TRUE, exist_ptime_col = "Step", exist_path_col = "Group"))`.
- `scmp.ob`: a fully processed `ScMaSigPro` object used by vignettes.
- `multi.lin.sce`: tradeSeq multifurcating simulation; three paths plus `root` in `group`, pseudotime in `Monocle3_Pseudotime`.

### Conventions

- Naming: pipeline steps are `sc.<verb>` (dots, mirroring maSigPro); accessors and plots are camelCase (`showSol`, `plotTrend`); constructors and helpers are snake_case (`as_scmp`, `pb_counts`); internals carry `@keywords internal`.
- Every exported function opens with `assertthat::assert_that(is(scmpObj, "ScMaSigPro"), msg = ...)` followed by column-existence checks with explicit messages.
- Code is styler / tidyverse formatted with 2-space indent. The `.Rproj` says 4 spaces, but the code and the CI `/style` command use styler defaults.
- Plots return ggplot objects; path colours come from `scmp_colors()` in R/utils.R; dplyr code uses `.data$col` / `!!sym()` to stay R CMD check clean.
- `parallel = TRUE` uses `parallel::mclapply` (fork), so it silently runs sequentially on Windows.

### Known quirks (verified in source, not yet fixed)

- `sc.t.fit()` overwrites its `p_value` argument with `scmpObj@Parameters@p_value` (R/sc.t.fit.R:85), so the argument has no effect; set `p_value` in `sc.p.vector()` instead.
- `plotDiagnostics()` composes its four panels with patchwork operators (`+`, `/`), but patchwork is not declared anywhere; it only works when patchwork is attached.
- `plotTrend()` uses `x^5` for the degree-6 term (R/plotTrend.R:152); harmless for `poly_degree <= 5`.
- The unexported `ScMaSigPro()` constructor in R/scmpGenerics.R passes `ParameterConfig =` for the `Parameters` slot. All real code uses `new("ScMaSigPro", ...)` directly.
- The `useWeights` blocks in `sc.p.vector` / `sc.t.fit` are hard-coded `FALSE` and reference a nonexistent `bin_size_colname` slot; dead code.
