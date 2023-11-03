#' Homogenize Bins of a Dataframe
#'
#' Adjust the bins of a dataframe based on a strategy that checks if the bin size is too small or too large
#' relative to the mean and standard deviation of the bin sizes. If a bin is too small, it will be merged
#' with the previous bin. If it's too large, it will be split into two equal parts.
#'
#' @param df A dataframe containing the bin data.
#' @param l_bound_col The name of the lower bound column as a string.
#' @param u_bound_col The name of the upper bound column as a string.
#' @param bin_size_col The name of the bin size column as a string.
#' @param binPTime_col The name of the column with Pseudotime
#'
#' @importFrom e1071 moment
#' @return A dataframe with adjusted bins.
#' @keywords internal
#

homogenize.bins <- function(df, l_bound_col, u_bound_col, bin_size_col, binPTime_col, verbose = TRUE) {
optimize_bins <- function(df) {
    mean_size <- mean(df[[bin_size_col]])
    sd_size <- sd(df[[bin_size_col]])

    if (verbose) {
        message(paste("Mean size:", round(mean_size), "SD size:", round(sd_size)))
    }

    # Initialize an empty dataframe with the same column names and types as df
    new_df <- df[0, ]

    for (i in 1:nrow(df)) {
        current_bin_size <- df[[bin_size_col]][i]

        if (current_bin_size < mean_size - sd_size && i != nrow(df)) {
            # Merge with the next bin because it is too small
            next_bin_size <- df[[bin_size_col]][i + 1]
            merged_bin_size <- current_bin_size + next_bin_size
            
            if(verbose){
                message(paste("Merging bin", i))
            }
            
            # Ensure the merged bin size is a whole number by adjusting the next bin size
            df[[bin_size_col]][i + 1] <- merged_bin_size
            df[[l_bound_col]][i + 1] <- df[[l_bound_col]][i]
        } else if (current_bin_size > mean_size + sd_size) {
            # Split the bin into smaller bins, ensuring whole numbers
            num_splits <- ceiling(current_bin_size / (mean_size + sd_size))
            split_size <- floor(current_bin_size / num_splits)
            remaining <- current_bin_size - split_size * num_splits
            
            if(verbose){
                message(paste("Splitting bin", i))
            }

            split_pts <- seq(df[[l_bound_col]][i], df[[u_bound_col]][i], length.out = num_splits + 1)
            for (j in 1:num_splits) {
                adjusted_split_size <- split_size
                if (j == num_splits) { # Add the remaining to the last split
                    adjusted_split_size <- adjusted_split_size + remaining
                }
                new_row <- setNames(data.frame(split_pts[j], split_pts[j + 1], adjusted_split_size, df[[binPTime_col]][i]),
                                    c(l_bound_col, u_bound_col, bin_size_col, binPTime_col))
                new_df <- rbind(new_df, new_row)
            }
        } else {
            # If the bin size is within acceptable range, keep it as it is
            new_df <- rbind(new_df, df[i, ])
        }
    }

    return(new_df)
}

  # Optimize bins iteratively until there's no change
  iteration_count <- 0
  previous_df <- df
  while (TRUE) {
    iteration_count <- iteration_count + 1
    if (verbose) {
      message(paste("Optimizing bin sizes, running iteration", iteration_count))
    }
    optimized_df <- optimize_bins(previous_df)
    if (identical(optimized_df, previous_df)) {
      if (verbose) {
        message("Optimization completed.")
      }
      break
    } else {
      previous_df <- optimized_df
    }
  }
  
  # Set rownames
  optimized_df[[binPTime_col]] <- seq(1:nrow(optimized_df)) 

  return(optimized_df)
}
