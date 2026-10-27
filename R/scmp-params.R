#' Arguments shared by several functions
#'
#' The arguments that several functions of the package take, documented
#' once and inherited by the help page of each function.
#'
#' @param scmpObj A [ScMaSigPro-class] object.
#' @param verbose Whether to print progress messages.
#' @param BPPARAM A `BiocParallelParam` object that runs the per-gene fits,
#'   such as `BiocParallel::MulticoreParam(2)` (Linux and macOS) or
#'   `BiocParallel::SnowParam(2)` (any system). The default,
#'   `BiocParallel::SerialParam()`, fits the genes one after another; create
#'   it with `progressbar = TRUE` for a progress bar. The results do not
#'   depend on the backend.
#' @param ptime_col The name of the column of the cell metadata that holds
#'   the pseudotime of each cell.
#' @param path_col The name of the column of the cell or bin metadata that
#'   holds the path label of each cell or bin.
#' @param bin_ptime_col The name of the column of the bin metadata that holds
#'   the binned pseudotime of each bin.
#' @param bin_col The name of the column of the bin metadata that holds the
#'   name of each bin.
#' @param bin_size_col The name of the column of the bin metadata that holds
#'   the number of cells of each bin.
#' @param bin_mem_col The name of the column of the bin metadata that holds
#'   the cells of each bin, separated by `|`.
#' @param assay_name The name of the assay of the `Sparse` slot that holds
#'   the counts of the cells.
#' @param aggregate How the counts of the cells of a bin are combined: `"sum"`
#'   or `"mean"` (rounded to whole numbers).
#' @param includeInflu Whether to keep the genes with influential
#'   observations (see [showInflu()]).
#' @param feature_id The name of a gene.
#'
#' @seealso [scMaSigPro-package]
#'
#' @name scmp-params
#' @keywords internal
NULL
