#' @title Rename the pseudotime and path columns of a SingleCellExperiment
#'
#' @description
#' `annotate_sce()` renames the columns `exist_ptime_col` and `exist_path_col`
#' of the cell metadata of `sce` to `ptime_col` and `path_col` when
#' `labels_exist = TRUE`, and returns `sce` unchanged otherwise.
#'
#' @param sce A `SingleCellExperiment`.
#' @param labels_exist Whether to rename the columns.
#' @param exist_ptime_col The column with the pseudotime.
#' @param exist_path_col The column with the path labels.
#' @inheritParams as_scmp
#'
#' @importFrom SingleCellExperiment colData
#'
#' @return `sce` with the renamed columns.
#'
#' @seealso [as_scmp()], which calls it.
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @keywords internal
annotate_sce <- function(sce,
                         ptime_col = "Pseudotime",
                         path_prefix = "Path",
                         root_label = "root",
                         path_col = path_prefix,
                         exist_ptime_col = NULL,
                         exist_path_col = NULL,
                         labels_exist = FALSE,
                         verbose = TRUE) {
    # Overwite the columns
    if (labels_exist) {
        .check(
            all(!is.null(exist_ptime_col) & !is.null(exist_path_col)),
            msg = paste(
                "Requested to set 'path_col' as", path_col,
                "with 'labels_exist' as TRUE. Please supply",
                "'exist_ptime_col' and 'exist_path_col'"
            )
        )

        if (verbose) {
            message(
                "Overwriting columns in cell level metadata, '", exist_path_col,
                "' is replaced by '", path_col, "' and '", exist_ptime_col,
                "' is replaced by '", ptime_col, "'."
            )
        }

        # Extract the cell metadata
        cell.meta <- as.data.frame(colData(sce))

        # Check columns
        .check(
            all(exist_ptime_col %in% colnames(cell.meta)),
            msg = paste0("'", exist_ptime_col, "' doesn't exist in colData.")
        )
        # Check columns
        .check(
            all(exist_path_col %in% colnames(cell.meta)),
            msg = paste0("'", exist_path_col, "' doesn't exist in colData.")
        )

        # Override
        names(cell.meta)[names(cell.meta) == exist_ptime_col] <- ptime_col

        # Overwite the columns
        names(cell.meta)[names(cell.meta) == exist_path_col] <- path_col

        # Update cell dataset with the updated cell metadata
        colData(sce) <- DataFrame(cell.meta)

        # Return
        return(sce)
    } else {
        if (verbose) {
            message("Skipping overwriting of colnames in cell level metadata")
        }
        # Return
        return(sce)
    }
}
