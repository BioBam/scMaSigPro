#' Simulated SingleCellExperiment Object
#'
#' A small simulated SingleCellExperiment Object created using Splatter.
#' This dataset contains 200 cells and 100 genes and is simulated to have
#' a bifurcating topology of the trajectory, useful for testing and development
#' in `scMaSigPro`. The dataset is stored as an `sce` object from the class
#' `SingleCellExperiment`
#'
#' @details
#' The `splat.sce` object was created using the `splatSimulatePaths` function
#' from the Splatter package. The following code was used for the simulation:
#' \preformatted{
#' # Load Required Packages
#' suppressPackageStartupMessages(library(splatter))
#' suppressPackageStartupMessages(library(scran))
#' suppressPackageStartupMessages(library(scuttle))
#' suppressPackageStartupMessages(library(scater))
#' suppressPackageStartupMessages(library(SingleCellExperiment))
#'
#' set.seed(123)
#'
#' # Simulate
#' splat.sim <- splatSimulatePaths(
#'    params = newSplatParams(
#'        batchCells = 200, nGenes = 100),
#'    group.prob = c(0.5, 0.5),
#'    path.nSteps = c(100, 100),
#'    de.prob = 0.3, de.facLoc = 0.2,
#'    path.from = c(0, 0), # Bifurcation
#'    verbose = FALSE)
#'
#' # Normalize
#'  splat.sim <- logNormCounts(splat.sim, assay.type = "counts")
#'
#' # Reduce Dimensions
#'  splat.sim <- runPCA(splat.sim, exprs_values = "logcounts", ncomponents = 2)
#'
#' # Visulize Steps and Groups
#'  plotPCA(splat.sim, colour_by = "Step")
#'  plotPCA(splat.sim, colour_by = "Group")
#'
#' # Create SCE and transfer data
#' sce <- SingleCellExperiment(list(counts = splat.sim@@assays@@data@@listData$counts))
#' sce@@colData <- splat.sim@@colData
#' rowData(sce) <- rowData(splat.sim)
#' reducedDims(sce) <- reducedDims(splat.sim)
#' splat.sim <- sce
#'
#' # Save
#'  save(splat.sim, file = "data/splat.sim.RData")
#'
#' # Compress
#' tools::resaveRdaFiles(paths = "data/")
#' }
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
