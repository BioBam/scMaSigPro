#' @title Wrapper for scMaSigPro Binning
#' @description
#' `squeeze()` performs a series of operations including data compression, creating pseudo-bulk-design, and creating pseudo-bulk counts.
#' It then generates a SingleCellExperiment object with updated slots.
#'
#' @param scmpObject An object.
#' @param pseudotime_colname Name of the column in `cell.metadata` generated using
#' \code{\link[SingleCellExperiment]{colData}} storing information for Pseudotime.
#' (Default is "Pseudotime")
#' @param path_colname Name of the column in `cell.metadata` generated using
#' \code{\link[SingleCellExperiment]{colData}} storing information for Path. 
#' (Default is `path_prefix`)
#' @param bin_method A character string (default = "Sturges"). The method to be
#' used to estimate the optimal number of bins.
#' @param drop.fac A numeric value (default = 0.5). The factor by which to
#' decrease the number of bins if the initial binning results in too many bins.
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' @param cluster_count_by A character string (default = "sum"). The method to
#' use to aggregate counts within each cluster.
#' @param assay_name The name of the assay (default: "counts").
#' @param binning A character string (deafult = "universal"). When set to
#' "individual", the bins are calculated per path iteratively.
#' @param additional_params Additional Parameters in a list
#' The name of the column in the design file that contains the binned data.
#'
#' @return Returns the object of 'scMaSigProClass' with updated slots.
#'
#' #' @examples
#' \dontrun{
#'
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
squeeze <- function(scmpObject,
                    pseudotime_colname = "Pseudotime",
                    path_colname = "Path",
                    bin_method = "Sturges",
                    drop.fac = 0.5,
                    verbose = TRUE,
                    cluster_count_by = "sum",
                    assay_name = "counts",
                    binning = "universal",
                    additional_params = NULL,
                    bin_pseudotime_colname = "scmp_binned_pseudotime") {
  # Extract the Cell Metadata
  cell.metadata <- as.data.frame(colData(scmpObject@sce))

  if (!is.null(additional_params)) {
    assert_that(is.list(additional_params),
      msg = "Please provide 'additional_params' as a named list. See details for more information"
    )

    assert_that(names(additional_params) %in% c("bin_pseudotime_colname", "bin_colname", "bin_size_colname", "bin_members_colname"),
      msg = "Incorrect additional parameters supplied , please see details."
    )
  }

  # Validate the data
  # pseudotime_colname
  assert_that((pseudotime_colname %in% colnames(cell.metadata)),
    msg = paste0("'", pseudotime_colname, "' ", "doesn't exit in cell.metadata.")
  )
  # path_colname
  assert_that((path_colname %in% colnames(cell.metadata)),
    msg = paste0("'", path_colname, "' ", "doesn't exit in cell.metadata.")
  )
  # drop.fac
  assert_that(drop.fac >= 0.3,
    msg = "Invalid value for 'drop.fac'. It should be between 0.3 and 1."
  )
  # Binning
  assert_that(all(binning %in% c("universal", "individual")),
    msg = "Allowed options for binning are 'Universal' and 'Individual'"
  )
  # Binning methods
  assert_that(
    all(
      bin_method %in% c("Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane", "Scott.Normal")
    ),
    msg = "Available binning methods are 'Freedman.Diaconis', 'Sqrt', 'Sturges', 'Rice', 'Doane', and 'Scott.Normal'"
  )
  # Count slot
  assert_that(
    all(
      assay_name %in% names(scmpObject@sce@assays@data@listData)
    ),
    msg = paste0("'", assay_name, "' ", "doesn't exit in scmpObject")
  )

  # Create Compression file
  tryCatch(
    expr = {
      compression.file <- entropy_discretize(
        cell_metadata = cell.metadata,
        pseudotime_colname = pseudotime_colname,
        path_colname = path_colname,
        bin_method = bin_method,
        drop.fac = drop.fac,
        verbose = verbose,
        binning = binning,
        bin_pseudotime_colname = bin_pseudotime_colname
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
        path_colname = path_colname,
        bin_pseudotime_colname = bin_pseudotime_colname
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
        counts = as.matrix(scmpObject@sce@assays@data@listData[[assay_name]]),
        bin_colname = "scmp_bin",
        pseudo_bulk_profile = compressed.cell.metadata,
        cluster_count_by = cluster_count_by
      )
    },
    error = function(e) {
      message(paste("Error message: ", e$message))
      stop("Unable to create pseudo-bulk-counts")
    }
  )

  # Create  a sce-object
  compressed.sce <- SingleCellExperiment(assays = list(bulk.counts = as(as.matrix(compressed.counts), "dgCMatrix")))
  colData(compressed.sce) <- DataFrame(compressed.cell.metadata)

  # Update Slot
  scmpObject@compress.sce <- compressed.sce

  # Return S3
  return(scmpObject)
}
