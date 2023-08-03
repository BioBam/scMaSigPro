#' Extract Slingshot components
#'
#' This function extracts the sling curves, sling time, and sling dimensions from a given SingleCellExperiment object
#' processed by the Slingshot package. The function also allows choosing the dimensionality reduction method to use.
#'
#' @param sling.sce A SingleCellExperiment object processed by Slingshot.
#' @param reduction_method A string indicating the dimensionality reduction method to use. Default is "umap".
#' @param verbose A logical value indicating whether to print detailed messages. Default is TRUE.
#'
#' @return A list containing the extracted sling curves, sling time, and sling dimensions, each as a dataframe.
#'
#' @examples
#' # Assuming 'sce' is a SingleCellExperiment object processed by Slingshot
#' result <- extract_slingshot_components(sce)
#' @export
#'
extract_slingshot_components <- function(sling.sce, reduction_method = "umap", verbose = TRUE) {
  # Convert the reduction method to upper case
  reduction_method <- toupper(reduction_method)

  # Extract Sling Curves
  sling.curves.df <- slingCurves(sling.sce, as.df = T)

  # Extract Sling Time
  sling.time.df <- slingPseudotime(sling.sce)

  # Extract Sling Dimensions
  sling.dim.df <- reducedDims(sling.sce)[[reduction_method]] %>%
    as.data.frame() %>%
    rename_with(~ c("dim_1", "dim_2"))

  return(list(
    sling.curves = sling.curves.df,
    sling.time = sling.time.df,
    sling.dim = sling.dim.df
  ))
}
