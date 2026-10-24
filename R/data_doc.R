#' Simulated bifurcating trajectory
#'
#' A `SingleCellExperiment` with the counts of 100 genes in 200 cells,
#' simulated with splatter (Zappia et al. 2017) along two paths that start
#' from the same point.
#'
#' @format A `SingleCellExperiment` with 100 genes (rows) and 200 cells
#' (columns):
#' \describe{
#'   \item{`counts` assay}{The simulated counts.}
#'   \item{`colData`}{`Cell`, the cell name; `Batch`, the batch (one,
#'     `"Batch1"`); `Group`, the path of the cell, a factor with the levels
#'     `"Path1"` and `"Path2"`; `ExpLibSize`, its expected library size;
#'     `Step`, how far along its path the cell is (1 to 100), used as the
#'     pseudotime; `sizeFactor`, the size factor of
#'     `scuttle::logNormCounts()`.}
#'   \item{`rowData`}{`Gene`, the gene name; `BaseGeneMean`, the base mean
#'     of the gene; `OutlierFactor`, its outlier factor; `GeneMean`, its
#'     mean after the outlier factor; `LocalDEFacPath1` and
#'     `LocalDEFacPath2`, its differential expression factor along each path;
#'     `DEFacPath1` and `DEFacPath2`, its factor at the end of each path;
#'     `SigmaFacPath1` and `SigmaFacPath2`, the factor of its non-linear
#'     change along each path (0 for a linear change).}
#'   \item{`reducedDims`}{`PCA`, the first two principal components of the
#'     log-normalised counts.}
#' }
#'
#' @source Simulated with `splatter::splatSimulatePaths()` (splatter 1.26.0,
#' `set.seed(123)`, 200 cells, 100 genes, two paths of 100 steps from the
#' origin, `de.prob = 0.3`, `de.facLoc = 0.2`), normalised with
#' `scuttle::logNormCounts()` and reduced with `scater::runPCA()`. The script
#' is `system.file("scripts", "make-splat-sim.R", package = "scMaSigPro")`.
#'
#' @references
#' Zappia L, Phipson B, Oshlack A (2017). Splatter: simulation of
#' single-cell RNA sequencing data. Genome Biology 18, 174.
#' \doi{10.1186/s13059-017-1305-0}
#'
#' @family datasets
#'
#' @usage data("splat.sim", package = "scMaSigPro")
#'
#' @examples
#' data("splat.sim", package = "scMaSigPro")
#' table(SummarizedExperiment::colData(splat.sim)$Group)
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
"splat.sim"

#' A complete analysis of splat.sim
#'
#' A [ScMaSigPro-class] object with every step of the workflow run on
#' [splat.sim], for the examples and the vignettes.
#'
#' @format A [ScMaSigPro-class] object:
#' \describe{
#'   \item{`Sparse`}{The 200 cells and 100 genes of [splat.sim], with the
#'     pseudotime in `Pseudotime` (`Step` of [splat.sim]) and the path in
#'     `Path` (`Group`).}
#'   \item{`Dense`}{The 16 bins of [sc.squeeze()] with its defaults, 8 per
#'     path, and their pseudo-bulk counts.}
#'   \item{`Design`}{The design of [sc.set.poly()] with
#'     `poly_degree = 1` and the offsets of [sc.p.vector()].}
#'   \item{`Profile`}{The test of the full model of every gene by
#'     [sc.p.vector()]: negative binomial family with theta 10, log link and
#'     size-factor offsets.}
#'   \item{`Estimate`}{The models selected by [sc.t.fit()] (backward
#'     selection).}
#'   \item{`Significant`}{The gene lists `Path1` and `Path2vsPath1` of
#'     [sc.filter()] with `rsq = 0.7` and `vars = "groups"`, and four
#'     clusters of their union from [sc.cluster.trend()].}
#'   \item{`Parameters`}{The settings of these steps.}
#' }
#'
#' @source Made by `system.file("scripts", "make-scmp-ob.R", package =
#' "scMaSigPro")` from [splat.sim].
#'
#' @references
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' @family datasets
#'
#' @usage data("scmp.ob", package = "scMaSigPro")
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' scmp.ob
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
"scmp.ob"

#' Simulated multifurcating trajectory
#'
#' A `SingleCellExperiment` with the counts of 501 genes in 1977 cells of a
#' simulated trajectory with three end points and a shared root, with the
#' pseudotime and the path labels of each cell.
#'
#' @format A `SingleCellExperiment` with 501 genes (rows) and 1977 cells
#' (columns):
#' \describe{
#'   \item{`counts` assay}{The simulated counts.}
#'   \item{`colData`}{`cell_id`, the cell name; `group`, the path label:
#'     `"root"` for the cells on all three lineages, `"path_Y_15"`,
#'     `"path_Y_18"` or `"path_Y_52"` for the cells on one, and
#'     `"path_Y_18|path_Y_15"` for the cells on these two;
#'     `Monocle3_Pseudotime`, the pseudotime of Monocle 3.}
#'   \item{`reducedDims`}{`TSNE`, the two dimensions of the t-SNE embedding.}
#' }
#' The `rowData` has no columns.
#'
#' @source The multifurcating simulation of Van den Berge et al. (2020),
#' made with dyntoy (Saelens et al. 2019): `multifurcating_4.rds` of
#' <https://github.com/statOmics/tradeSeqPaper>. Monocle 3 (Cao et al. 2019)
#' normalised the counts, embedded the cells with PCA and t-SNE, learnt the
#' principal graph and gave the pseudotime from the root node `Y_51`. Each
#' cell is on the lineages whose shortest path from the root, on that graph,
#' holds its nearest node, as in the Monocle 3 steps of the tradeSeq
#' vignette; three lineages are kept and the cells on none of them are
#' left out. The script is
#' `system.file("scripts", "make-multi-lin-sce.R", package = "scMaSigPro")`.
#'
#' @references
#' Van den Berge K, Roux de Bézieux H, Street K, Saelens W, Cannoodt R,
#' Saeys Y, Dudoit S, Clement L (2020). Trajectory-based differential
#' expression analysis for single-cell sequencing data. Nature
#' Communications 11, 1201. \doi{10.1038/s41467-020-14766-3}
#'
#' Saelens W, Cannoodt R, Todorov H, Saeys Y (2019). A comparison of
#' single-cell trajectory inference methods. Nature Biotechnology 37(5),
#' 547-554. \doi{10.1038/s41587-019-0071-9}
#'
#' Cao J, Spielmann M, Qiu X, Huang X, Ibrahim DM, Hill AJ, Zhang F,
#' Mundlos S, Christiansen L, Steemers FJ, Trapnell C, Shendure J (2019).
#' The single-cell transcriptional landscape of mammalian organogenesis.
#' Nature 566(7745), 496-502. \doi{10.1038/s41586-019-0969-x}
#'
#' @family datasets
#'
#' @usage data("multi.lin.sce", package = "scMaSigPro")
#'
#' @examples
#' data("multi.lin.sce", package = "scMaSigPro")
#' table(SummarizedExperiment::colData(multi.lin.sce)$group)
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
"multi.lin.sce"
