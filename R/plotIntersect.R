#' @title Generate UpSet Plot on Intersection of Significant Genes from scmpObject
#'
#' @param scmpObj An object of class ScMaSigPro
#' @param min_intersection_size minimal number of observations in an intersection
#' for it to be included.
#' @param width_ratio ratio of the overall set size width to intersection matrix
#' width.
#' @param keep_empty_groups whether empty sets should be kept (including sets
#' which are only empty after filtering by size)
#' @param show_sets_size the overall set sizes plot, e.g. from upset_set_size()
#' (FALSE to hide)
#'
#' @return An UpSet plot visualizing the intersections of significant genes across pathways.
#' @importFrom S4Vectors isEmpty
#' @importFrom ComplexUpset upset intersection_matrix intersection_size upset_set_size
#' @importFrom RColorConesa colorConesa
#'
#' @export
#'
plotIntersect <- function(scmpObj, min_intersection_size = 2,
                          keep_empty_groups = TRUE,
                          width_ratio = 0.1, show_sets_size = FALSE) {
  # Check the data
  assert_that(
    is(scmpObj, "ScMaSigPro"),
    msg = "Please supply an object of the class 'ScMaSigPro'"
  )

  # Check if siggenes results exist for groups
  assert_that(!isEmpty(scmpObj@sig.genes@sig.genes),
    msg = "'sig.genes@Summary' slot is empty, please run 'sc.get.siggenes'"
  )

  gene_list <- scmpObj@sig.genes@sig.genes
  # Create a unique list of all genes
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
    show_sets_size <- upset_set_size()
  }

  # Create Upset
  p <- upset(
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

      intersection_matrix(
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
      "Intersection size" = intersection_size(
        counts = TRUE,
        mapping = aes(fill = "bars_color")
      )
      + scale_fill_manual(values = c("bars_color" = col_pal[2]), guide = "none")
    )
  ) + ggtitle("Intersection of features among paths") +
    theme(legend.position = "none")

  # return plot
  return(p)

  # Perform UpSet plot
  # upset(
  #     fromList(gene.list),
  #     main.bar.color = "#F58A53",
  #     matrix.color = "#15918A",
  #     line.size = 1.5,
  #     point.size = 3,
  #     shade.color = "purple",
  #     text.scale = 1.5,
  #     sets.x.label = "Number of Features",
  #     sets.bar.color = "#EE446F"
  # )
}
