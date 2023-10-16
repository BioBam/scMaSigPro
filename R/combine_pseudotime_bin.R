#' Combine pseudotime and bin data
#'
#' This function takes a pseudotime data frame (path.frame) and a bin data frame (bin_table) and combines them.
#' It subsets path.frame based on the bin_table 'from' and 'to' values and adds this bin_table information to path.frame.
#' The function is parallelized to speed up the process.
#'
#' @param path.frame A data frame containing pseudotime data. It must include a 'Step' column.
#' @param bin_table A data frame containing bin data. It must include 'from', 'to', 'bin_size', and 'binnedTime' columns.
#' @param time.col A data frame containing bin data. It must include 'from', 'to', 'bin_size', and 'binnedTime' columns.
#'
#' @return A data frame which is a combined version of path.frame and bin_table.
#'
#' @importFrom parallel detectCores mclapply
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'path.frame' and 'bin_table' are available and have the correct format
#' combined_data <- combine_pseudotime_bin(path.frame, bin_table)
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @keywords internal
combine_pseudotime_bin <- function(path.frame, bin_table, time.col) {
  # Get the number of cores
  num_cores <- detectCores()

  # Use mclapply to apply an anonymous function in parallel
  processed_design_table <- mclapply(seq_len(nrow(bin_table)), function(i, time_col = time.col) {
    # Subset path.frame where Step is within the current bin_table from and to
    subset_df <- path.frame[path.frame[[time_col]] >= bin_table$from[i] & path.frame[[time_col]] <= bin_table$to[i], ]

    # If the subset is not empty, add bin_table information
    if (nrow(subset_df) > 0) {
      subset_df$from <- bin_table$from[i]
      subset_df$to <- bin_table$to[i]
      subset_df$bin_size <- bin_table$bin_size[i]
      subset_df$binnedTime <- bin_table$binnedTime[i]

      return(subset_df)
    }

    return(data.frame()) # Return an empty data.frame if no rows match
  }, mc.cores = num_cores)


  # Combine all list elements into a single data frame
  processed_design_table <- do.call(rbind, processed_design_table)

  return(processed_design_table)
}
