#' @title Create a ScMaSigPro object from a SingleCellExperiment
#'
#' @description
#' `as_scmp()` creates a [ScMaSigPro-class] object from a
#' `SingleCellExperiment` whose cell metadata (`colData`) holds a pseudotime
#' column and a path column. To start from a count matrix and a cell table,
#' use [create_scmp()]; [assign_paths()] and [paths_from_graph()] build the
#' two columns from the output of a trajectory tool.
#'
#' @param object A `SingleCellExperiment` with the counts of the cells.
#' @param from The class of `object`; only `"sce"` is supported.
#' @param path_prefix The default of `path_col`, kept in the parameters.
#' @param root_label The label of the root cells, kept in the parameters.
#' @param additional_params A list: with `labels_exist = TRUE`, the columns
#'   `exist_ptime_col` and `exist_path_col` of the cell metadata are renamed
#'   to `ptime_col` and `path_col`.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @return A [ScMaSigPro-class] object with the cells in the `Sparse` slot.
#'
#' @examples
#' data("splat.sim", package = "scMaSigPro")
#' scmp <- as_scmp(splat.sim,
#'     additional_params = list(
#'         labels_exist = TRUE, exist_ptime_col = "Step",
#'         exist_path_col = "Group"
#'     ),
#'     verbose = FALSE
#' )
#' scmp
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @family input helpers
#' @inheritParams scmp-params
#'
#' @export
as_scmp <- function(object, from = "sce",
                    path_prefix = "Path",
                    root_label = "root",
                    ptime_col = "Pseudotime",
                    path_col = path_prefix,
                    verbose = TRUE,
                    additional_params = list(
                        labels_exist = FALSE
                    )) {
    # Check Conversion Type
    .check(identical(from, "sce"),
        msg = paste(
            "'from' must be \"sce\": as_scmp() converts SingleCellExperiment",
            "objects only. Monocle3 'cell_data_set' objects (from = \"cds\")",
            "are no longer converted; pass the counts and a cell table with",
            "pseudotime and path columns to create_scmp(), as shown in the",
            "'Common Root Cells' vignette."
        )
    )

    # Validate S4
    .check(is(object, "SingleCellExperiment"),
        msg = paste(
            "'object' must be a SingleCellExperiment. To start from a count",
            "matrix and a table of cell metadata, use create_scmp()."
        )
    )
    # Names and labels, stored in the Parameters slot
    .check_strings(
        path_prefix = path_prefix, root_label = root_label,
        ptime_col = ptime_col, path_col = path_col
    )
    # An object saved by an older SingleCellExperiment is updated here, as
    # colData<- does with labels_exist = TRUE, so both routes store the same
    # object.
    object <- updateObject(object)

    # Check and validate additional parameters
    if (!is.null(additional_params)) {
        .check(is.list(additional_params),
            msg = paste(
                "Please provide 'additional_params' as a named list.",
                "See details for more information"
            )
        )
        # Check additional parameters
        .check(
            all(names(additional_params) %in% c(
                "exist_ptime_col", "exist_path_col", "labels_exist"
            )),
            msg = paste(
                "Allowed additional parameters for sce are 'exist_ptime_col',",
                "'exist_path_col', and 'labels_exist'."
            )
        )
    } else {
        additional_params <- list(
            labels_exist = NULL,
            exist_ptime_col = NULL,
            exist_path_col = NULL
        )
    }

    # Without labels_exist, the columns keep their names
    labels_exist <- additional_params[["labels_exist"]]
    if (is.null(labels_exist)) {
        labels_exist <- FALSE
    }

    if (verbose) {
        message("Supplied object: SingleCellExperiment object")
    }
    # Annotate the sce object
    annotated_sce <- annotate_sce(
        sce = object,
        ptime_col = ptime_col,
        path_prefix = path_prefix,
        root_label = root_label,
        path_col = path_col,
        exist_ptime_col = additional_params[["exist_ptime_col"]],
        exist_path_col = additional_params[["exist_path_col"]],
        labels_exist = labels_exist,
        verbose = verbose
    )

    # Create Object
    scmpObj <- new("ScMaSigPro",
        Sparse = annotated_sce,
        Dense = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
    )

    # Update the Parametersaram slot
    scmpObj@Parameters@ptime_col <- ptime_col
    scmpObj@Parameters@root_label <- root_label
    scmpObj@Parameters@path_prefix <- path_prefix
    scmpObj@Parameters@path_col <- path_col
    validObject(scmpObj)

    # Return Object
    return(scmpObj)
}
