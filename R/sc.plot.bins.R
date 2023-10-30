#' sc.plot.bins: Plot Bin Sizes Across Binned Time and Paths
#'
#' This function generates plots to visualize the compressed data from a SingleCellExperiment object.
#' It produces a bar plot and, optionally, a tile (heatmap) plot to display the bin sizes across different
#' binned time intervals and paths.
#'
#' @param scmpObj A SingleCellExperiment object with an additional slot '@compress.sce' that contains compression information.
#' @param add_tile A logical flag indicating whether to add a tile (heatmap) plot alongside the bar plot. Default is TRUE.
#' @param bar.title A title of the barplot
#' @param tile.title A title of the heatmap
#'
#' @return A bar plot and, optionally, a tile (heatmap) plot, visualizing the bin sizes across different binned time and paths.
#' If add_tile is TRUE, returns a combined ggplot object with both plots; otherwise, only the bar plot is printed.
#'
#'
#' @export
# sc.plot.bins <- function(scmpObj, add_tile = T, bar.title = "Bin size across binned pseudotime",
#                          tile.title = "Bin size across paths along binned pseudotime") {
#   # Extract the compressed data
#   compression.info <- as.data.frame(colData(scmpObj@compress.sce))
#
#   conesa_colors <- getConesaColors()[c(T, F)][c(1:length(unique(compression.info[["path"]])))]
#   names(conesa_colors) <- unique(unique(compression.info[["path"]]))
#
#
#   bar <- ggplot(compression.info, aes(x = factor(binnedTime), y = bin.size)) +
#     geom_bar(stat = "identity", aes(fill = path), position = "dodge") +
#     geom_line(aes(group = path, color = path), position = position_dodge(0.9)) +
#     geom_point(aes(color = path), position = position_dodge(0.9)) +
#     ggtitle(bar.title) +
#     scale_fill_manual(values = conesa_colors) +
#     xlab("Binned Time") +
#     ylab("Bin Size") +
#     theme_minimal()
#
#   if (add_tile) {
#     tile <- ggplot(compression.info, aes(x = factor(binnedTime), y = path)) +
#       geom_tile(aes(fill = bin.size)) +
#       scale_fill_gradient(low = "#FDA3D1", high = "#FDC659") +
#       geom_text(aes(label = sprintf("%d", round(bin.size, 1))), vjust = 1) +
#       ggtitle(tile.title) +
#       xlab("Binned Time") +
#       ylab("Path") +
#       theme_minimal()
#
#
#     ggarrange(bar, tile, nrow = 2)
#   } else {
#     print(bar)
#   }
# }
