#' @title Create ScMaSigPro Object
#'
#' @description
#' `create_scmp()` initializes a ScMaSigPro object with the given counts,
#' cell level metadata, and other optional parameters.
#'
#' @param counts A matrix containing the raw expression counts.
#' @param cell_data A data frame containing the cell level metadata.
#' @param feature_data Optional data frame of gene metadata, one row per gene
#' of `counts` (named like them, or with automatic row names). It becomes the
#' `rowData()` of the `Sparse` slot, and of the `Dense` slot when this
#' function builds it.
#' @param use_as_bin A logical indicating to run MaSigPro analysis without
#' pseudotime binning. Defaults to FALSE.
#' @param bin_counts Optional pseudo-bulk counts of bins made outside
#' scMaSigPro: a matrix with the genes of `counts` as rows, in the same
#' order, and one column per bin. Requires `bin_cell_data`; not with
#' `use_as_bin = TRUE`.
#' @param bin_cell_data A data frame with one row per bin, named like the
#' columns of `bin_counts`, holding the binned pseudotime in `ptime_col` and
#' the path in `path_col`. A `scmp_bin_size` column gives the number of cells
#' of each bin; without it the sizes are `NA`.
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values. (Default is "Pseudotime")
#' @param path_col A character string representing the column name for branching
#' path assignment. (Default is `path_prefix`)
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with aligned pseudotime.
#'
#' @details
#' The `Dense` slot holds the bins. By default it is empty and
#' [sc.squeeze()] bins the cells. With `use_as_bin = TRUE` every cell is a
#' bin. With `bin_counts` and `bin_cell_data` it holds bins made outside
#' scMaSigPro, and the workflow continues with [sc.set.poly()]; such bins
#' have no member cells, so [pb_counts()] and [sc.restruct()] do not apply.
#'
#' @importFrom S4Vectors DataFrame
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @export

# Create ScMaSigPro
create_scmp <- function(counts,
                        cell_data,
                        feature_data = NULL,
                        bin_counts = NULL,
                        bin_cell_data = NULL,
                        ptime_col,
                        path_col,
                        use_as_bin = FALSE) {
    # Validation Checks
    .check(ncol(counts) == nrow(cell_data),
        msg = paste(
            "Number of cells in raw-counts and cell-level-metadata are",
            "different."
        )
    )

    .check(all(colnames(counts) == rownames(cell_data)),
        msg = paste(
            "Rownames of raw-counts and cell-level-metadata are different."
        )
    )

    use_bins <- !is.null(bin_counts) || !is.null(bin_cell_data)
    if (use_bins) {
        .check(!is.null(bin_counts) && !is.null(bin_cell_data),
            msg = "Please provide both 'bin_counts' and 'bin_cell_data'."
        )
        .check(!use_as_bin,
            msg = "'bin_counts' cannot be combined with 'use_as_bin = TRUE'."
        )
        .check(ncol(bin_counts) == nrow(bin_cell_data),
            msg = paste(
                "Number of bins in bin_counts and bin_cell_data are different."
            )
        )
        .check(identical(colnames(bin_counts), rownames(bin_cell_data)),
            msg = paste(
                "Colnames of bin_counts and rownames of bin_cell_data are",
                "different."
            )
        )
        .check(identical(rownames(bin_counts), rownames(counts)),
            msg = "The rows of bin_counts must be the genes of counts."
        )
    }
    if (!is.null(feature_data)) {
        feature_data <- as.data.frame(feature_data)
        .check(nrow(feature_data) == nrow(counts),
            msg = "feature_data must have one row per gene of counts."
        )
        if (.row_names_info(feature_data) < 0) {
            # Automatic row names: the rows are the genes of counts
            rownames(feature_data) <- rownames(counts)
        }
        .check(identical(rownames(feature_data), rownames(counts)),
            msg = "The rownames of feature_data must be the genes of counts."
        )
        feature_data <- DataFrame(feature_data)
    }

    .check_strings(path_col = path_col, ptime_col = ptime_col)

    # Validate if path_col is present in cell_data
    .check(path_col %in% colnames(cell_data),
        msg = paste0("Column name '", path_col, "' not found in cell_data.")
    )

    # Validate if ptime_col is present in cell_data
    .check(ptime_col %in% colnames(cell_data),
        msg = paste0("Column name '", ptime_col, "' not found in cell_data.")
    )
    if (use_bins) {
        for (col in c(path_col, ptime_col)) {
            .check(col %in% colnames(bin_cell_data),
                msg = paste0(
                    "Column name '", col, "' not found in bin_cell_data."
                )
            )
        }
    }

    # Create Single-Cell Experiment Object
    sparse_tmp <- SingleCellExperiment(
        list(counts = counts),
        colData = DataFrame(cell_data), rowData = feature_data
    )

    # Initate scMaSigPro
    scmpObj <- new("ScMaSigPro",
        Sparse = sparse_tmp,
        Dense = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
    )
    rm(sparse_tmp)

    # Use as bin
    if (use_as_bin) {
        # Set bin size coulmn name
        cell_data[["scmp_bin_size"]] <- as.numeric(1)

        # Create sparse
        sparse_tmp <- SingleCellExperiment(
            list(bulk.counts = counts),
            colData = DataFrame(cell_data), rowData = feature_data
        )

        # Transfer Data
        scmpObj@Dense <- sparse_tmp
        rm(sparse_tmp)

        # Update the slots
        scmpObj@Parameters@bin_ptime_col <- ptime_col
    }

    # Bins made outside scMaSigPro: the binned pseudotime is in ptime_col
    if (use_bins) {
        bin_cell_data <- as.data.frame(bin_cell_data)
        bin_size_col <- scmpObj@Parameters@bin_size_col
        if (!bin_size_col %in% colnames(bin_cell_data)) {
            bin_cell_data[[bin_size_col]] <- NA_real_
        }
        if (!inherits(bin_counts, "dgCMatrix")) {
            bin_counts <- as(as.matrix(bin_counts), "dgCMatrix")
        }
        scmpObj@Dense <- SingleCellExperiment(
            list(bulk.counts = bin_counts),
            colData = DataFrame(bin_cell_data), rowData = feature_data
        )
        scmpObj@Parameters@bin_ptime_col <- ptime_col
    }

    # Update the slots
    scmpObj@Parameters@ptime_col <- ptime_col
    scmpObj@Parameters@path_col <- path_col
    validObject(scmpObj)

    return(scmpObj)
}
