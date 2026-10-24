#' @title Plot the sizes of the bins
#'
#' @description
#' `plotBinTile()` draws a tile for each bin, by path and binned pseudotime,
#' coloured and labelled by its number of cells.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 .data aes geom_text geom_tile ggplot labs
#'   scale_fill_gradient theme_minimal xlab ylab
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @family plots
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' plotBinTile(scmp.ob)
#'
#' @export
plotBinTile <- function(scmpObj,
                        path_col = scmpObj@Parameters@path_col,
                        bin_size_col = scmpObj@Parameters@bin_size_col,
                        bin_ptime_col = scmpObj@Parameters@bin_ptime_col) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check whether the compression data exist or not
    compression.info <- as.data.frame(colData(scmpObj@Dense))

    # Check if values are binned
    .check(nrow(as.data.frame(colData(scmpObj@Dense))) >= 1,
        msg = "No binning information found. Please run 'sc.squeeze()', first."
    )

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
