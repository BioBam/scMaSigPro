#' @title Generate UpSet Plot
#'
#' @description
#' Generate UpSet Plot on Intersection of Significant Genes from scMaSigPro
#' object. It is a wrapper `UpSetR::upset`.
#'
#' @importFrom S4Vectors isEmpty
#' @importFrom RColorConesa colorConesa
#' @importFrom utils packageVersion
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param min_intersection_size Minimal number of observations in an intersection
#' for it to be included.
#' @param width_ratio Ratio of the overall set size width to intersection matrix
#' width.
#' @param keep_empty_groups Whether empty sets should be kept (including sets
#' which are only empty after filtering by size)
#' @param show_sets_size The overall set sizes plot, e.g. from upset_set_size()
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return upset object for 'UpSetR'.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
plotIntersect <- function(scmpObj,
                          min_intersection_size = 2,
                          keep_empty_groups = TRUE,
                          width_ratio = 0.1,
                          show_sets_size = FALSE,
                          verbose = TRUE) {
  # Check the data
  assertthat::assert_that(
    is(scmpObj, "ScMaSigPro"),
    msg = "Please supply an object of the class 'ScMaSigPro'"
  )

  # Check if siggenes results exist for groups
  assertthat::assert_that(!isEmpty(scmpObj@Significant@genes),
    msg = "'sig.genes@Summary' slot is empty, please run 'sc.get.siggenes'"
  )

  # Check if package is installed
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(paste0("Package '", package, "' is not installed. Please install it first."))
  } else {
    if (verbose) {
      message(paste0("Using '", package, "' for UpSet plot."))
    }
  }

  gene_list <- scmpObj@Significant@genes

  # Create list object
  upset_r_gene_list <- UpSetR::fromList(gene_list)

  # Create Plot
  p <- UpSetR::upset(
    upset_r_gene_list,
    main.bar.color = "#F58A53",
    matrix.color = "#15918A",
    line.size = 1.5,
    cutoff = min_intersection_size,
    empty.intersections = keep_empty_groups,
    point.size = 3,
    shade.color = "purple",
    text.scale = 1.5,
    sets.x.label = "Number of Features",
    sets.bar.color = "#EE446F"
  )
  return(p)
}
