#' @description
#' scMaSigPro finds genes whose expression changes along pseudotime and
#' differs between the branching paths of a single-cell trajectory. The
#' cells of each path are pseudo-bulked into pseudotime bins, and the
#' polynomial generalised linear model of maSigPro is fitted to the binned
#' counts, by default with a negative binomial family and size-factor
#' offsets (Srivastava et al. 2024).
#'
#' @section Workflow:
#' 1. [as_scmp()] or [create_scmp()]: a [ScMaSigPro-class] object from the
#'    counts and one pseudotime value and one path label per cell.
#' 2. [sc.squeeze()]: pseudotime bins and pseudo-bulk counts for each path;
#'    [sc.restruct()] adds the bins of a shared root to every path.
#' 3. [sc.set.poly()]: the polynomial design.
#' 4. [sc.p.vector()]: the test of each gene's full model, which selects the
#'    genes with a non-flat profile.
#' 5. [sc.t.fit()]: stepwise selection of the significant terms of each of
#'    these genes.
#' 6. [sc.filter()]: lists of significant genes; [sc.cluster.trend()]
#'    clusters them.
#' 7. [plotTrend()], [plotTrendCluster()], [plotIntersect()]: plots;
#'    [scmpResults()]: the results as a table with one row per gene.
#'
#' @section Input:
#' The workflow needs one pseudotime value and one path label per cell.
#' [assign_paths()] and [paths_from_graph()] build them from the lineage
#' weights or the graph of a trajectory tool; the "Input contract" section
#' of [assign_paths()] describes the format.
#'
#' @section Documentation:
#' The vignettes and these help pages are also on the package website,
#' <https://spriyansh.github.io/scMaSigPro/>:
#' * Quick Start Guide:
#'   <https://spriyansh.github.io/scMaSigPro/articles/Basic-Workflow.html>
#' * Common Root Cells:
#'   <https://spriyansh.github.io/scMaSigPro/articles/Common-Root-Cells.html>
#' * scMaSigProClass:
#'   <https://spriyansh.github.io/scMaSigPro/articles/scMaSigPro-Class.html>
#' * Bridging maSigPro Analysis:
#'   <https://spriyansh.github.io/scMaSigPro/articles/scMaSigPro-maSigPro.html>
#'
#' @references
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' Conesa A, Nueda MJ, Ferrer A, Talón M (2006). maSigPro: a method to
#' identify significantly differential expression profiles in time-course
#' microarray experiments. Bioinformatics 22(9), 1096-1102.
#' \doi{10.1093/bioinformatics/btl056}
#'
#' Nueda MJ, Tarazona S, Conesa A (2014). Next maSigPro: updating maSigPro
#' bioconductor package for RNA-seq time series. Bioinformatics 30(18),
#' 2598-2602. \doi{10.1093/bioinformatics/btu333}
#'
#' @keywords internal
"_PACKAGE"
