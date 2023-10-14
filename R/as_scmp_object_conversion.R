#' Convert cell_data_set or SingleCellExperiment to scMaSigProClass
#'
#' This function converts a cell_data_set object from Monocle or a
#' SingleCellExperiment object from slingshot to an instance of the
#' scMaSigProClass.
#'
#' @param object An S4 object of class 'cell_data_set' or 'SingleCellExperiment'.
#' @param from Character string specifying the class of 'object'. Can be either
#' "cds" for "cell_data_set" or "sce" for "SingleCellExperiment".
#' @param path_prefix Prefix used to annoate the paths, default is "Path".
#' @param root_label Label used to describe root cells,  default is "root".
#' @param pseudotime_colname Name of the column in cell.metadata storing
#' information for Pseudotime.
#' @param path_colname Name of the column in cell.metadata storing information
#' about annoated paths.
#' @param additional_params A list of additional parameters, see details.
#' 
#' @details Additional parameters 
#' 
#' 
#' @return An instance of the 'scMaSigProClass'.
#' @examples
#' \dontrun{
#' scmpObj <- as_scmp(object,
#'   from = "cds", path_prefix = "Path",
#'   root_label = "root", path_colname = path_prefix,
#'   verbose = TRUE
#' )
#' }
#'
#' @importFrom monocle3 cell_data_set
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
    msg = ("Currently, accepted objects are from 'cell_data_set' and 'SingleCellExperiment' classes")
  )

  # Validate S4
  assert_that(
    all(isS4(object) & all(is(object, "cell_data_set") | is(object, "SingleCellExperiment"))),
    msg = "Please provide object from one of the class 'Monocle3/cell_data_set', 'SingleCellExperiment/SCE'"
  )

  # Check and validate additional parameters
  if (!is.null(additional_params)) {
    assert_that(is.list(additional_params),
      msg = "Please provide 'additional_params' as a named list. See details for more information"
    )
    # Check additional parameters
    if (from == "cds") {
      assert_that(names(additional_params) %in% c("reduction_method"),
        msg = "Additional Parameters for CDS are 'reduction_method'."
      )
    } else if (from == "sce") {
      assert_that(all(names(additional_params) %in% c("existing_pseudotime_colname", "existing_path_colname", "overwrite_labels")),
        msg = "Additional Parameters for SCE are 'existing_pseudotime_colname',
                  'existing_path_colname', 'overwrite_labels'."
      )
    }
  }else{
      if (from == "cds") {
          additional_params <- list(reduction_method = "umap")
      } else if (from == "sce") {
          additional_params <- list(existing_pseudotime_colname = NULL,
                                    existing_path_colname = NULL)
      }
  }

  # FLow control
  if (is(object)[1] == "SingleCellExperiment") {
    if (verbose) {
      message("Supplied object: SingleCellExperiment object")
    }
      # Annotate the monocel3 Object
      annotated_sce <- annotate_sce(sce = object,
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
    
    # Return Object
    return(scmpObj)
    
    
  } else if (is(object)[1] == "cell_data_set") {
    if (verbose) {
      message("Supplied object: cell_data_set object from Monocle3")
    }
    # Annotate the monocel3 Object
    annotated_cds <- annotate_monocle3_cds(cds,
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
    return(scmpObj)
  }
}
