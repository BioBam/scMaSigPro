#' @title Convert a 'SingleCellExperiment' object to a ScMaSigPro object.
#'
#' @description
#' `as_scmp()` converts a `SingleCellExperiment` object, whose cell metadata
#' (`colData`) holds a pseudotime column and a branching-path column, to a
#' \code{\link{ScMaSigPro}} object. To start from a count matrix and a table of
#' cell metadata, for example with pseudotime and path labels extracted from a
#' trajectory tool such as Monocle3 or slingshot, use [create_scmp()]; the
#' "Common Root Cells" vignette shows the steps for a Monocle3 `cell_data_set`.
#' [assign_paths()] and [paths_from_graph()] build the two columns from
#' lineage weights or from a trajectory graph.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @param object An S4 object of class `SingleCellExperiment`.
#' @param from Character string specifying the class of 'object'. Only "sce"
#' (`SingleCellExperiment`) is supported; Monocle3 `cell_data_set` objects
#' ("cds") are no longer converted. (Default is "sce").
#' @param path_prefix Prefix used to annotate the paths. (Default is "Path").
#' @param root_label Label used to annotate root cells. (Default is "root").
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values. (Default is "Pseudotime")
#' @param path_col A character string representing the column name for branching
#' path assignment. (Default is `path_prefix`)
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' @param additional_params A named list of additional parameters:
#' `labels_exist`, `exist_ptime_col` and `exist_path_col`. With
#' `labels_exist = TRUE`, the existing `colData` columns `exist_ptime_col` and
#' `exist_path_col` are renamed to `ptime_col` and `path_col`. See examples.
#'
#' @return An object of class \code{\link{ScMaSigPro}}.
#'
#' @examples
#' # Load Library
#' # library(scMaSigPro)
#' # Step-1: Load a dataset for testing
#' # This dataset is available as part of the package
#' # It is simulated with splatter
#' data("splat.sim", package = "scMaSigPro")
#'
#' # Step-2: Convert to ScMaSigPro Object
#' # Here, we convert the sce object to an scMaSigPro object
#' scmp.sce <- as_scmp(
#'     object = splat.sim, from = "sce",
#'     verbose = FALSE,
#'     additional_params = list(
#'         labels_exist = TRUE,
#'         exist_ptime_col = "Step",
#'         exist_path_col = "Group"
#'     )
#' )
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
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
