#' @title Plot Bin Sizes Across Binned Time and Paths
#'
#' @description
#' This function generates plots to visualize the Dense slot cell metadata
#' from a ScMaSigPro object. It produces tile plot to display the bin sizes
#' across different binned time intervals and paths.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param path_col A character string representing the column name for branching
#' path assignment in 'Sparse' or 'Dense' slot.
#' @param bin_size_col A character string representing the name of the column in
#' which bin sizes per bin are stored. (Default is "scmp_bin_size").
#' @param bin_ptime_col A character string representing the column name
#' for binned Pseudotime values in 'Dense' data.
#' (Default is "scmp_binned_pseudotime").
#'
#' @return ggplot2 plot object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
plotBinTile <- function(scmpObj,
                        path_col = scmpObj@Parameters@path_col,
                        bin_size_col = scmpObj@Parameters@bin_size_col,
                        bin_ptime_col = scmpObj@Parameters@bin_ptime_col) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Check whether the compression data exist or not
  compression.info <- as.data.frame(colData(scmpObj@Dense))

  # Check for extended data
  if (nrow(compression.info) < 1) {
    compression.info <- as.data.frame(colData(scmpObj@Sparse))
  }

  # Check if values are binned
  assert_that(nrow(compression.info) >= 1,
    msg = "Please run 'sc.squeeze()' first."
  )

  # get conesa colors
  conesa_colors <- getConesaColors()[c(TRUE, FALSE)][c(1:length(unique(compression.info[[path_col]])))]
  names(conesa_colors) <- unique(unique(compression.info[[path_col]]))

  # Create plot data
  plt.data <- data.frame(
    pTime = as.factor(compression.info[[bin_ptime_col]]),
    bin = compression.info[[bin_ptime_col]],
    path = compression.info[[path_col]],
    binSize = compression.info[[bin_size_col]]
  )

  # Create plot
  tile <- ggplot(plt.data, aes(x = .data$pTime, y = .data$path)) +
    geom_tile(aes(fill = .data$binSize)) +
    scale_fill_gradient(low = "#FDA3D1", high = "#FDC659") +
    geom_text(aes(label = sprintf("%d", round(.data$binSize, 1))), vjust = 1) +
    labs(fill = bin_size_col) +
    xlab(bin_ptime_col) +
    ylab(path_col) +
    theme_minimal()

  return(tile)
}
