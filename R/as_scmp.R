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
#' @param align_pseudotime Whether to automatically align two different pseudotimes.
#' See \code{\link{align.pseudotime}} for more details. (Default is FALSE).
#' @param interactive Whether to use the shiny application to select paths. (Default is TRUE).
#' @param annotation_colname Column name in cell level metadata containing cell type
#' annotations. (Default is "cell_type").
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
#' # Load Library
#' # library(scMaSigPro)
#' # Step-1: Load a dataset for testing
#' # This dataset is available as part of the package
#' # It is simulated with splatter
#' data("Sim2Path", package = "scMaSigPro")
#'
#' # Step-2: Convert to ScMaSigPro Object
#' # Here, we convert the SCE object to an scMaSigPro object
#' scmp.sce <- as_scmp(
#'   object = sim.sce, from = "sce",
#'   align_pseudotime = TRUE,
#'   verbose = FALSE,
#'   additional_params = list(
#'     labels_exist = TRUE,
#'     existing_pseudotime_colname = "Step",
#'     existing_path_colname = "Group"
#'   )
#' )
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
                    align_pseudotime = FALSE,
                    path_colname = path_prefix,
                    annotation_colname = "cell_type",
                    verbose = TRUE,
                    interactive = TRUE,
                    additional_params = list(
                      labels_exist = FALSE
                    )) {
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
      assert_that(names(additional_params) %in% c("reduction_method", "labels_exist", "align_pseudotime_method"),
        msg = "Allowed additional parameters for 'cds' (cds/CellDataSet) are
        'reduction_method', and 'labels_exist','align_pseudotime_method'."
      )
    } else if (from == "sce") {
      assert_that(all(names(additional_params) %in% c("existing_pseudotime_colname", "existing_path_colname", "labels_exist")),
        msg = "Allowed additional parameters for SCE are 'existing_pseudotime_colname',
                  'existing_path_colname', and 'labels_exist'."
      )
    }
  } else {
    # If additional_params == NULL for 'cds'
    if (from == "cds") {
      additional_params <- list(reduction_method = "umap")
    } else if (from == "sce") {
      additional_params <- list(
        labels_exist = NULL,
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
      labels_exist = additional_params[["labels_exist"]],
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
    # annotated_cds <- annotate_monocle3_cds(object,
    #   reduction_method = additional_params[["reduction_method"]],
    #   path_prefix = path_prefix,
    #   root_label = root_label,
    #   path_colname = path_colname,
    #   pseudotime_colname = pseudotime_colname,
    #   verbose = verbose
    # )
    if (interactive) {
      scmpObj <- selectPath.m3(
        cdsObj = object,
        annotation_col = annotation_colname,
        pseudotime_col = pseudotime_colname,
        path_col = path_colname,
        redDim = additional_params[["reduction_method"]]
      )
    } else {
      stop("Only support for interactive for now")
    }
    # return(scmpObj)
    if (align_pseudotime) {
      scmpObj <- align.pseudotime(
        scmpObj = scmpObj,
        method = "rescale",
        pseudotime_col = pseudotime_colname,
        path_col = path_colname,
        verbose = verbose
      )
    }

    # Print

    return(scmpObj)
  }
}
