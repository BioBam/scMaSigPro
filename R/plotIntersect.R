#' @title Generate UpSet Plot
#'
#' @description
#' Generate UpSet Plot on Intersection of Significant Genes from scMaSigPro
#' object. It is a wrapper around `ComplexUpset::upset` and `UpSetR::upset`.
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
#' @param package Which package to use for the UpsetPlot. Options are 'ComplexUpset'
#' or 'UpSetR' (Default).
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return ggplot2 plot object for 'ComplexUpset' or upset object for 'UpSetR'.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
plotIntersect <- function(scmpObj,
                          package = "UpSetR",
                          min_intersection_size = 2,
                          keep_empty_groups = TRUE,
                          width_ratio = 0.1,
                          show_sets_size = FALSE,
                          verbose = TRUE) {
  # Check the data
  assert_that(
    is(scmpObj, "ScMaSigPro"),
    msg = "Please supply an object of the class 'ScMaSigPro'"
  )

  # Check if siggenes results exist for groups
  assert_that(!isEmpty(scmpObj@Significant@genes),
    msg = "'sig.genes@Summary' slot is empty, please run 'sc.get.siggenes'"
  )

  # Check for possible options
  assert_that(package %in% c("ComplexUpset", "UpSetR"),
    msg = "Please provide a valid package name for UpSet plot. Options are 'ComplexUpset' or 'UpSetR'"
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

  if (package == "UpSetR") {
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
  } else {
    # Check version of the ggplot2
    if (packageVersion("ggplot2") >= "3.5.0") {
      warning("Please downgrade the ggplot2 to '>= 3.5.0' to use 'ComplesUpset'. We will support the latest version in future.
          Visit:'https://github.com/krassowski/complex-upset/issues/195' for more details.")
    } else {
      # # Create a unique list of all genes
      all_genes <- unique(unlist(gene_list))

      # Initialize the data frame
      gene_df <- data.frame(gene = all_genes)

      # Add columns for each pathway
      for (pathway in names(gene_list)) {
        gene_df[[pathway]] <- gene_df$gene %in% gene_list[[pathway]]
      }

      # Binarize Variables and set factors
      gene_df[, -1] <- lapply(gene_df[, -1], function(x) as.integer(x))

      # Get conesa colours
      col_pal <- colorConesa(3)

      if (show_sets_size) {
        show_sets_size <- ComplexUpset::upset_set_size()
      }

      # Create Upset
      p <- ComplexUpset::upset(
        data = gene_df,
        intersect = colnames(gene_df)[-1],
        width_ratio = width_ratio,
        min_size = min_intersection_size,
        keep_empty_groups = keep_empty_groups,
        name = "Vars",
        # wrap=FALSE,
        set_sizes = show_sets_size,
        # stripes=c('deepskyblue1'),
        matrix = (

          ComplexUpset::intersection_matrix(
            geom = geom_point(
              shape = "square",
              size = 3.5
            ),
            segment = geom_segment(
              linetype = "dotted",
              color = col_pal[1]
            )
          )
          + scale_color_manual(
              values = c("TRUE" = col_pal[1], "FALSE" = col_pal[3]),
              # labels=c('TRUE'='yes', 'FALSE'='no'),
              breaks = c("TRUE", "FALSE")
            )
        ),
        base_annotations = list(
          "Intersection size" = ComplexUpset::intersection_size(
            counts = TRUE,
            mapping = aes(fill = "bars_color")
          )
          + scale_fill_manual(values = c("bars_color" = col_pal[2]), guide = "none")
        )
      ) + ggtitle("Intersection of features among paths") +
        theme(legend.position = "none", legend.title = element_text(hjust = 0.5))

      # return plot
      return(p)
    }
  }
}
