#' @title Convert 'Cell Dataset' or 'SingleCellExperiment' object to scmpClass
#'
#' @description
#' `as_scmp()` converts a cds/CellDataSet object from Monocle3 or a SingleCellExperiment
#' object from Slingshot to an instance of the scmpClass object.
#'
#' @param object An S4 object of class `cds/CellDataSet` or `SingleCellExperiment`.
#' @param from Character string specifying the class of 'object'. Use "cds" for
#' `cds/CellDataSet` class and "sce" for `SingleCellExperiment` class.
#' @param path_prefix Prefix used to annotate the paths. (Default is "Path").
#' @param root_label Label used to annotate root cells. (Default is "root").
#' @param pseudotime_colname Name of the column in `cell.metadata` storing
#' Pseudotime values. It is generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is "Pseudotime").
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path. It is generated using `colData` from the \pkg{SingleCellExperiment} package.
#' (Default is `path_prefix`)
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' @param additional_params A named list of additional parameters. See details.
#'
#' @details Additional Details
#'
#' @return An instance of the 'scmpClass'.
#'
#' @seealso `colData` from the \pkg{SingleCellExperiment} package, `new_cell_data_set`
#' function in \pkg{monocle3}
#'
#' @examples
#' \dontrun{
#' scmpObj <- as_scmp(object,
#'   from = "cds", path_prefix = "Path",
#'   root_label = "root",
#'   path_colname = path_prefix,
#'   verbose = TRUE
#' )
#' }
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom assertthat assert_that
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
as_scmp <- function(object, from = "cds",
                    path_prefix = "Path",
                    root_label = "root",
                    pseudotime_colname = "Pseudotime",
                    path_colname = path_prefix,
                    verbose = TRUE,
                    additional_params = list(overwrite_labels = TRUE)) {
  # Check Conversion Type
  assert_that(from %in% c("cds", "sce"),
    msg = ("Currently, accepted options in the 'from' parameter are 'cds'
    ('cds/CellDataSet' object) and 'sce' ('SingleCellExperiment').")
  )

  # Validate S4
  assert_that(
    all(isS4(object) & all(is(object, "cell_data_set") | is(object, "SingleCellExperiment"))),
    msg = "Please provide object from one of the class 'cds/CellDataSet',
    or 'SingleCellExperiment/SCE'."
  )

  # Check and validate additional parameters
  if (!is.null(additional_params)) {
    assert_that(is.list(additional_params),
      msg = "Please provide 'additional_params' as a named list.
      See details for more information"
    )
    # Check additional parameters
    if (from == "cds") {
      assert_that(names(additional_params) %in% c("reduction_method", "overwrite_labels"),
        msg = "Allowed additional parameters for 'cds' (cds/CellDataSet) are
        'reduction_method', and 'overwrite_labels'."
      )
    } else if (from == "sce") {
      assert_that(all(names(additional_params) %in% c("existing_pseudotime_colname", "existing_path_colname", "overwrite_labels")),
        msg = "Allowed additional parameters for SCE are 'existing_pseudotime_colname',
                  'existing_path_colname', and 'overwrite_labels'."
      )
    }
  } else {
    # If additional_params == NULL for 'cds'
    if (from == "cds") {
      additional_params <- list(reduction_method = "umap")
    } else if (from == "sce") {
      additional_params <- list(
        existing_pseudotime_colname = NULL,
        existing_path_colname = NULL
      )
    }
  }

  # if 'sce' is selected
  if (is(object)[1] == "SingleCellExperiment") {
    if (verbose) {
      message("Supplied object: SingleCellExperiment object")
    }
    # Annotate the sce object
    annotated_sce <- annotate_sce(
      sce = object,
      pseudotime_colname = pseudotime_colname,
      path_prefix = path_prefix,
      root_label = root_label,
      path_colname = path_colname,
      existing_pseudotime_colname = additional_params[["existing_pseudotime_colname"]],
      existing_path_colname = additional_params[["existing_path_colname"]],
      overwrite_labels = additional_params[["overwrite_labels"]],
      verbose = verbose
    )

    # Create Object
    scmpObj <- new("scMaSigProClass",
      sce = annotated_sce,
      compress.sce = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
    )

    # Update the AddParams slot
    scmpObj@addParams@pseudotime_colname <- pseudotime_colname
    scmpObj@addParams@root_label <- root_label
    scmpObj@addParams@path_prefix <- path_prefix
    scmpObj@addParams@path_colname <- path_colname

    if (verbose) {
      print(scmpObj)
    }

    # Return Object
    return(scmpObj)
  } else if (is(object)[1] == "cell_data_set") {
    if (verbose) {
      message("Supplied object: cell_data_set object from Monocle3")
    }

    # Add extraction
    if (is.null(additional_params[["reduction_method"]])) {
      additional_params[["reduction_method"]] <- "umap"
    }

    # Annotate the monocel3 Object
    annotated_cds <- annotate_monocle3_cds(object,
      reduction_method = additional_params[["reduction_method"]],
      path_prefix = path_prefix,
      root_label = root_label,
      path_colname = path_colname,
      pseudotime_colname = pseudotime_colname,
      verbose = verbose
    )

    # Convert to sce abd then to scMaSigPro Class
    scmpObj <- new("scMaSigProClass",
      sce = annotated_cds,
      compress.sce = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
    )

    # Update Slots
    scmpObj@addParams@pseudotime_colname <- pseudotime_colname
    scmpObj@addParams@root_label <- root_label
    scmpObj@addParams@path_prefix <- path_prefix
    scmpObj@addParams@path_colname <- path_colname

    # Print
    if (verbose) {
      print(scmpObj)
    }
    return(scmpObj)
  }
}
