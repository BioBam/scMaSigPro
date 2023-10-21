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
                    pseudotime_colname = scmpObject@addParams@pseudotime_colname,
                    path_colname = scmpObject@addParams@path_colname,
                    bin_method = "Sturges",
                    drop.fac = 0.5,
                    verbose = TRUE,
                    cluster_count_by = "sum",
                    assay_name = "counts",
                    binning = "universal",
                    additional_params = NULL,
                    bin_pseudotime_colname = scmpObject@addParams@bin_pseudotime_colname) {
  if (!is.null(additional_params)) {
    assert_that(is.list(additional_params),
      msg = "Please provide 'additional_params' as a named list. See details for more information"
    )

    assert_that(names(additional_params) %in% c("bin_pseudotime_colname", "bin_colname", "bin_size_colname", "bin_members_colname"),
      msg = "Incorrect additional parameters supplied , please see details."
    )
  }

  # Create Compression file
  tryCatch(
    expr = {
      scmpObject <- entropy_discretize(
        scmpObject = scmpObject,
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
      scmpObject <- make.pseudobulk.design(
        scmpObject = scmpObject,
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
      scmpObject <- make.pseudobulk.counts(
        scmpObject = scmpObject, # as.matrix(scmpObject@sce@assays@data@listData[[assay_name]]),
        bin_colname = "scmp_bin",
        cluster_count_by = cluster_count_by
      )
    },
    error = function(e) {
      message(paste("Error message: ", e$message))
      stop("Unable to create pseudo-bulk-counts")
    }
  )

  # Return S3
  return(scmpObject)
}
