#' @title Convert 'Cell Dataset' or 'SingleCellExperiment' object to ScMaSigPro
#' object.
#'
#' @description
#' `as_scmp()` converts a cds/CellDataSet object from Monocle3 or a
#' SingleCellExperiment #' object to an instance of the scmpClass object.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom assertthat assert_that
#'
#' @param object An S4 object of class `cds/CellDataSet` or `SingleCellExperiment`.
#' @param from Character string specifying the class of 'object'. Use "cds" for
#' `cds/CellDataSet` class and "sce" for `SingleCellExperiment` class.
#' @param path_prefix Prefix used to annotate the paths. (Default is "Path").
#' @param root_label Label used to annotate root cells. (Default is "root").
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values. (Default is "Pseudotime")
#' @param path_col A character string representing the column name for branching
#' path assignment. (Default is `path_prefix`)
#' @param align_pseudotime Whether to automatically align two different pseudotimes.
#' See \code{\link{align_pseudotime}} for more details. (Default is FALSE).
#' @param interactive Whether to use the shiny application to select branching
#' paths. (Default is TRUE).
#' @param anno_col A character string representing the column name for cell level
#' metadata containing cell level annotations. (Default is "cell_type").
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' @param additional_params A named list of additional parameters. See examples.
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
#'   object = splat.sim, from = "sce",
#'   align_pseudotime = TRUE,
#'   verbose = FALSE,
#'   additional_params = list(
#'     labels_exist = TRUE,
#'     exist_ptime_col = "Step",
#'     exist_path_col = "Group"
#'   )
#' )
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
as_scmp <- function(object, from = "cds",
                    path_prefix = "Path",
                    root_label = "root",
                    ptime_col = "Pseudotime",
                    align_pseudotime = FALSE,
                    path_col = path_prefix,
                    anno_col = "cell_type",
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
    or 'SingleCellExperiment/sce'."
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
      assert_that(all(names(additional_params) %in% c("exist_ptime_col", "exist_path_col", "labels_exist")),
        msg = "Allowed additional parameters for sce are 'exist_ptime_col',
                  'exist_path_col', and 'labels_exist'."
      )
    }
  } else {
    # If additional_params == NULL for 'cds'
    if (from == "cds") {
      additional_params <- list(reduction_method = "umap")
    } else if (from == "sce") {
      additional_params <- list(
        labels_exist = NULL,
        exist_ptime_col = NULL,
        exist_path_col = NULL
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
      ptime_col = ptime_col,
      path_prefix = path_prefix,
      root_label = root_label,
      path_col = path_col,
      exist_ptime_col = additional_params[["exist_ptime_col"]],
      exist_path_col = additional_params[["exist_path_col"]],
      labels_exist = additional_params[["labels_exist"]],
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

    if (interactive) {
      scmpObj <- m3_select_path(
        cds = object,
        anno_col = anno_col,
        ptime_col = ptime_col,
        path_col = path_col,
        latent_space = additional_params[["reduction_method"]]
      )
    } else {
      stop("Only support for interactive for now")
    }
    # return(scmpObj)
    if (align_pseudotime) {
      scmpObj <- align_pseudotime(
        scmpObj = scmpObj,
        method = "rescale",
        ptime_col = ptime_col,
        path_col = path_col,
        verbose = verbose
      )
    }

    # Print

    return(scmpObj)
  }
}
