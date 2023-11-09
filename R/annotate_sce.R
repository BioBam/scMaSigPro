#' @title Annotate 'SingleCellExperiment' class object with pseudotime and path information.
#'
#' @description
#' `annotate_sce()` annotates a SingleCellExperiment class object with pseudotime
#' and path information in its `cell.metadata` generated using `colData` from the
#' \pkg{SingleCellExperiment} package.
#'
#' @param sce A SingleCellExperiment object to be annotated.
#' @param pseudotime_colname Name of the column in `cell.metadata`  storing
#' information for Pseudotime. It is generated using `colData` from the
#' \pkg{SingleCellExperiment} package. (Default is "Pseudotime").
#' @param path_prefix Prefix used to annotate the paths. (Default is "Path").
#' @param root_label Label used to annotate root cells. (Default is "root").
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path. It is generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is `path_prefix`).
#' @param existing_pseudotime_colname The name of an existing pseudotime column
#' to be replaced (if not NULL).
#' @param existing_path_colname The name of an existing path column to be replaced
#' (if not NULL).
#' @param labels_exist Logical, should existing column names be overwritten
#' if they already exist? (default is TRUE).
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return A SingleCellExperiment object with updated cell metadata.
#'
#' @details Additional Details
#'
#' @seealso SingleCellExperiment class object, `colData` from the \pkg{SingleCellExperiment} package.
#'
#' @examples
#' # Annotate a SingleCellExperiment object with pseudotime and path information
#' \donttest{
#' annotated_sce <- annotate_sce(sce,
#'   pseudotime_colname = "Pseudotime",
#'   path_colname = "Path"
#' )
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom assertthat assert_that
#'
#' @export
annotate_sce <- function(sce,
                         pseudotime_colname = "Pseudotime",
                         path_prefix = "Path",
                         root_label = "root",
                         path_colname = path_prefix,
                         existing_pseudotime_colname = NULL,
                         existing_path_colname = NULL,
                         labels_exist = FALSE,
                         verbose = TRUE) {
  # Overwite the columns
  if (labels_exist) {
    assert_that(
      all(!is.null(existing_pseudotime_colname) & !is.null(existing_path_colname)),
      msg = paste("Requested to set 'path_colname' as", path_colname, "with 'labels_exist' as TRUE. Please supply,
                  'existing_pseudotime_colname' and 'existing_path_colname', through 'additional_params'")
    )

    if (verbose) {
      message(paste0(
        "Overwritting columns in cell.level.metadata, '", existing_path_colname, "' is replaced by '", path_colname,
        "' and '", existing_pseudotime_colname, "' is replaced by '", pseudotime_colname, "'."
      ))
    }

    # Extract the cell metadata
    cell.meta <- as.data.frame(colData(sce))

    # Check columns
    assert_that(
      all(existing_pseudotime_colname %in% colnames(cell.meta)),
      msg = paste("'", existing_pseudotime_colname, "', doesn't exist in cell.metadata.")
    )
    # Check columns
    assert_that(
      all(existing_path_colname %in% colnames(cell.meta)),
      msg = paste("'", existing_path_colname, "', doesn't exist in cell.metadata.")
    )

    # Override
    names(cell.meta)[names(cell.meta) == existing_pseudotime_colname] <- pseudotime_colname

    # Overwite the columns
    names(cell.meta)[names(cell.meta) == existing_path_colname] <- path_colname

    # Update cell dataset with the updated cell metadata
    sce@colData <- DataFrame(cell.meta)

    # Return
    return(sce)
  } else {
    if (verbose) {
      message(paste("Skipping overwritting of colnames in cell.level.metadata requested"))
    }
    # Return
    return(sce)
  }
}
