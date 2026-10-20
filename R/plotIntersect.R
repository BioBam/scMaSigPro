#' @title Generate UpSet Plot
#'
#' @description
#' Generate UpSet Plot on Intersection of Significant Genes from scMaSigPro
#' object. It is a wrapper `UpSetR::upset`, so drawing the plot needs the
#' suggested package 'UpSetR'; `return = TRUE` does not.
#'
#' @importFrom S4Vectors isEmpty
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param min_intersection_size Minimal number of observations in an intersection
#' for it to be included.
#' @param width_ratio Ratio of the overall set size width to intersection matrix
#' width.
#' @param keep_empty_groups Whether empty sets should be kept (including sets
#' which are only empty after filtering by size)
#' @param show_sets_size The overall set sizes plot, e.g. from upset_set_size()
#' @param return If set to true, it will return dataframe from the UpSetR::fromList().
#' (Default is TRUE)
#'
#' @return upset object for 'UpSetR'.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
plotIntersect <- function(scmpObj,
                          min_intersection_size = 1,
                          keep_empty_groups = FALSE,
                          width_ratio = 0.1,
                          show_sets_size = FALSE,
                          return = FALSE) {
    # Check the data
    .check_scmp(scmpObj)

    # Check if siggenes results exist for groups
    .check(!isEmpty(scmpObj@Significant@genes),
        msg = "'Significant' slot is empty, please run 'sc.filter()' first."
    )

    if (!keep_empty_groups) {
        keep_empty_groups <- NULL
    }

    gene_list <- scmpObj@Significant@genes

    # Create list object
    upset_r_gene_list <- fromListWithNames(gene_list)

    if (return) {
        return(as.data.frame(upset_r_gene_list))
    } else {
        .check(requireNamespace("UpSetR", quietly = TRUE),
            msg = paste(
                "plotIntersect() draws the plot with the 'UpSetR' package;",
                "install it with install.packages('UpSetR'), or use",
                "return = TRUE to get the data."
            )
        )
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
            sets.bar.color = "#EE446F",
            keep.order = TRUE,
            order.by = "freq"
        )
        return(p)
    }
}
