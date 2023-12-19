#' @title Plot Bin Sizes Across Binned Time and Paths
#'
#' @description
#' This function generates plots to visualize the dense cell metadata from a scmp
#' object. It produces tile plot to display the bin sizes across different
#' binned time intervals and paths.
#'
#' @param scmpObj A scmp class object with an additional slot '@dense' that
#' contains compression information.
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path.
#' @param bin_size_colname A title of the barplot
#' @param bin_pseudotime_colname description
#'
#' @return A tile plot made with `geom_tile()`, visualizing the bin sizes across
#' different binned time and paths.
#' @export
plotBinTile <- function(scmpObj,
                        path_colname = scmpObj@param@path_colname,
                        bin_size_colname = scmpObj@param@bin_size_colname,
                        bin_pseudotime_colname = scmpObj@param@bin_pseudotime_colname) {
  # Check Object Validity
  assert_that(is(scmpObj, "scmp"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Check whether the compression data exist or not
  compression.info <- as.data.frame(colData(scmpObj@dense))

  # Check for extended data
  if (nrow(compression.info) < 1) {
    compression.info <- as.data.frame(colData(scmpObj@sparse))
  }

  # Check if values are binned
  assert_that(nrow(compression.info) >= 1,
    msg = "Please run 'sc.squeeze()' first."
  )

  # get conesa colors
  conesa_colors <- getConesaColors()[c(TRUE, FALSE)][c(1:length(unique(compression.info[[path_colname]])))]
  names(conesa_colors) <- unique(unique(compression.info[[path_colname]]))

  # Create plot data
  plt.data <- data.frame(
    pTime = as.factor(compression.info[[bin_pseudotime_colname]]),
    bin = compression.info[[bin_pseudotime_colname]],
    path = compression.info[[path_colname]],
    binSize = compression.info[[bin_size_colname]]
  )

  # Create plot
  tile <- ggplot(plt.data, aes(x = .data$pTime, y = .data$path)) +
    geom_tile(aes(fill = .data$binSize)) +
    scale_fill_gradient(low = "#FDA3D1", high = "#FDC659") +
    geom_text(aes(label = sprintf("%d", round(.data$binSize, 1))), vjust = 1) +
    labs(fill = bin_size_colname) +
    xlab(bin_pseudotime_colname) +
    ylab(path_colname) +
    theme_minimal()

  return(tile)
}
