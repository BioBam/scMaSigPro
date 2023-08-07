#' The 'squeeze' function
#'
#' This function performs a series of operations including data compression, creating pseudo-bulk-design, and creating pseudo-bulk counts.
#' It then generates a SingleCellExperiment object with updated slots.
#'
#' @param scmp.ob An object.
#' @param time.col The name of the time column.
#' @param path.col The name of the path column.
#' @param method Method for data compression (default: "Sturges").
#' @param drop.fac The factor by which to drop data (default: 0.5).
#' @param verbose Logical. If TRUE, detailed messages are printed (default: TRUE).
#' @param cluster.count.by How to count clusters (default: "sum").
#' @param assay.name The name of the assay (default: "counts").
#'
#' @return Returns the object 'scmp.ob' with updated slots.
#' @export
#'
#' @examples
#' # Insert an example of how to use the function here.
#'
squeeze <- function(scmp.ob,
                    time.col,
                    path.col,
                    method = "Sturges",
                    drop.fac = 0.5,
                    verbose = TRUE,
                    cluster.count.by = "sum",
                    assay.name = "counts") {
  # Extract the Cell Metadata
  cell.metadata <- as.data.frame(colData(scmp.ob@sce))

  # Create Compression file
  tryCatch(
    expr = {
      compression.file <- entropy_discretize(cell.metadata,
        time_col = time.col,
        path_col = path.col,
        method = method,
        drop.fac = drop.fac,
        verbose = verbose
      )
    },
    error = function(e) {
      message(paste("Error message: ", e$message))
      stop("Unable to compress data")
    }
  )

  # Make New cell metaData
  tryCatch(
    expr = {
      compressed.cell.metadata <- make.pseudobulk.design(
        design.file = compression.file,
        pathCol = path.col,
        binnedCol = "binnedTime"
      )
    },
    error = function(e) {
      message(paste("Error message: ", e$message))
      stop("Unable to create pseudo-bulk-design")
    }
  )

  # Create Compressed Data
  tryCatch(
    expr = {
      compressed.counts <- make.pseudobulk.counts(
        counts = scmp.ob@sce@assays@data@listData[[assay.name]],
        cluster_member_col = "cluster.members",
        bin_col = "bin",
        pseudo_bulk_profile = compressed.cell.metadata,
        cluster.count.by = "sum"
      )
    },
    error = function(e) {
      message(paste("Error message: ", e$message))
      stop("Unable to create pseudo-bulk-design")
    }
  )

  # Create  a sce-object
  compressed.sce <- SingleCellExperiment(assays = list(bulk.counts = as(as.matrix(compressed.counts), "dgCMatrix")))
  colData(compressed.sce) <- DataFrame(compressed.cell.metadata)

  # Update Slot
  scmp.ob@compress.sce <- compressed.sce

  # Return S3
  return(scmp.ob)
}
