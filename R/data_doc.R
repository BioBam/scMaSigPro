#' @title Simulated SingleCellExperiment Object
#'
#' @description
#' A small simulated SingleCellExperiment Object created using Splatter.
#' This dataset contains 200 cells and 100 genes and is simulated to have
#' a bifurcating topology of the trajectory, useful for testing and development
#' in `scMaSigPro`. The dataset is stored as an `sce` object from the class
#' `SingleCellExperiment`
#'
#' @details
#' The `splat.sim` object was created with `splatter::splatSimulatePaths()`
#' (bifurcating topology, 200 cells, 100 genes), log-normalised with
#' `scuttle::logNormCounts()` and reduced with PCA. The complete script is
#' shipped with the package:
#' `system.file("scripts", "make-splat-sim.R", package = "scMaSigPro")`.
#'
#' This simulation creates a dataset with 100 genes and 200 cells, designed to mimic
#' a bifurcating trajectory typically observed in cellular differentiation.
#'
#' @usage
#' # Loading
#' data("splat.sim", package = "scMaSigPro")
#'
#' @format
#' An object of class `SingleCellExperiment` with 100 gene and 200 cells.
#'
#' @source
#' Simulated using the `Splatter` (1.26.0) package.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
"splat.sim"


#' @title scMaSigPro Object with results
#'
#' @description
#' A small scMaSigPro object created from simulated data. Please follow Vignette
#' `scMaSigPro: Quick Start Guide` for more details.
#'
#' @format An object of class `ScMaSigPro` built from `splat.sim` by the Quick
#' Start workflow (polynomial degree 1, negative binomial with log size-factor
#' offsets). The script that regenerates it is
#' `system.file("scripts", "make-scmp-ob.R", package = "scMaSigPro")`.
#'
#' @usage
#' # Loading
#' data("scmp.ob", package = "scMaSigPro")
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
"scmp.ob"

#' @title Simulated Multifurcating Trajectory SingleCellExperiment Object
#'
#' @description A simulated multifurcating trajectory derived from the tradeSeq
#' article and analysed with Monocle3, useful for testing and development of
#' trajectory analysis methods. The dataset is stored as a
#' `SingleCellExperiment` with raw `counts`, a t-SNE embedding in
#' `reducedDims(x)$TSNE`, and `colData` columns `cell_id`, `group` (hard branch
#' assignment: `root`, `path_Y_15`, `path_Y_18`, `path_Y_52` and the shared
#' segment `path_Y_18|path_Y_15`) and `Monocle3_Pseudotime`.
#'
#' @details
#'
#' The `multi.lin.sce` object was created using the following steps:
#'
#' 1. The multifurcating trajectory data was downloaded from the `tradeSeq` article.
#' 2. The raw counts, cell metadata, and gene metadata were extracted and transformed into a `Monocle3` `cell_data_set` object.
#' 3. The data was preprocessed using PCA for dimensionality reduction and normalized using a log transformation.
#' 4. Both t-SNE and UMAP were used for dimensionality reduction, with t-SNE embeddings stored in the UMAP slots to enable graph learning.
#' 5. Cells were clustered and a principal graph was learned on the data.
#' 6. Pseudotime was inferred based on the learned graph.
#'
#' The complete script is shipped with the package:
#' `system.file("scripts", "make-multi-lin-sce.R", package = "scMaSigPro")`.
#' This dataset includes expression data, cell metadata, and gene metadata, and it is structured to facilitate the application of various trajectory analysis methods.
#'
#' @format A `SingleCellExperiment` with 501 genes and 1977 cells.
#'
#' @name multi.lin.sce
#' @docType data
#' @keywords datasets
#'
#' @references
#' TradeSeq package: https://github.com/statOmics/tradeSeq
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
"multi.lin.sce"
