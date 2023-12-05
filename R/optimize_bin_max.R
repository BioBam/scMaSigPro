#' @title Homogenize Bins of a Dataframe
#'
#' @description
#' Adjust the bins of a dataframe based on a specified maximum bin size threshold.
#' It checks if the bin size is too large relative to this threshold.
#' If a bin is too large, it will be split into two equal parts until all bins are
#' below the maximum allowed size. Optionally, bins too small can be merged with
#' adjacent bins based on the provided method.
#'
#' @param bin_table A dataframe containing the bin data.
#' @param max_allowed The maximum allowed size for a bin.
#' @param verbose Logical; if TRUE, prints detailed output.
#' @param time_vector A vector representing time or a sequence that bins refer to.
#' @param lbound The name of the lower bound column in `bin_table`.
#' @param ubound The name of the upper bound column in `bin_table`.
#' @param bin The name of the bin identifier column in `bin_table`.
#' @param bin_size The name of the bin size column in `bin_table`.
#' @param method The method for handling small bins: 'merge' to merge with previous or next bin,
#' 'drop' to remove small bins, or 'ignore' to leave small bins as they are.
#' @param drop The threshold below which a bin is considered too small and subject to the method.
#'
#' @return A dataframe with adjusted bins. The structure of the dataframe will be the same as `bin_table`
#' with updated values for the bin size and bounds.
#' @keywords internal
optimize_bin_max <- function(bin_table, max_allowed, verbose = TRUE,
                             time_vector, lbound, ubound, bin, bin.size, method,
                             drop) {
  # Initiate an empty uniform bin
  uniform_bin_df <- data.frame(matrix(NA, ncol = ncol(bin_table)))
  colnames(uniform_bin_df) <- colnames(bin_table)

  # Binning one-by one
  for (i in c(1:nrow(bin_table))) {
    # histo_bin
    current_bin_df <- c(unlist(bin_table[i, , drop = FALSE]))

    # Get size of the bins
    current_bin_size <- current_bin_df[[bin.size]]

    # Check if the bin size is big
    if (current_bin_size > max_allowed) {
      # get pseudotime
      pTime.for.big.interval <- time_vector[(time_vector >= current_bin_df[[lbound]] & time_vector <= current_bin_df[[ubound]])]

      # Estimate how many splits are required
      potential_splits <- ceiling(estBinSize(
        time_vector = pTime.for.big.interval, nPoints = length(pTime.for.big.interval),
        drop_fac = drop, bin_method = method
      ))

      # Run the extraction
      small_splitted_bin_df <- extract_interval(
        time.vector = pTime.for.big.interval,
        nBins = potential_splits,
        bin = bin, bin.size = bin.size, lbound = lbound, ubound = ubound
      )
      # Validation
      if (verbose) {
        message(paste("Splitting bin with", current_bin_size, "cells into", nrow(small_splitted_bin_df), "with sizes as", paste(small_splitted_bin_df[[bin.size]], collapse = ", ")))
      }

      # Add new Row
      uniform_bin_df <- rbind(uniform_bin_df, small_splitted_bin_df)
    } else if (current_bin_size <= max_allowed) {
      uniform_bin_df <- rbind(uniform_bin_df, current_bin_df)
      if (verbose) {
        if (current_bin_size != 0) {
          message(paste("Skipping bin with size", current_bin_size))
        }
      }
    }
  }
  uniform_bin_df <- uniform_bin_df[-1, ]

  # Correct the rows
  rownames(uniform_bin_df) <- NULL

  # Return the new data
  return(uniform_bin_df)
}
