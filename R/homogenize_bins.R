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
#' @return A dataframe with adjusted bins.
#' @keywords internal
# 
homogenize_bins <- function(df, l_bound_col, u_bound_col, bin_size_col, binPTime_col, verbose = TRUE) {
    # Calculate mean and standard deviation of bin sizes
    mean_size <- mean(df[[bin_size_col]])
    sd_size <- sd(df[[bin_size_col]])
    
    optimize_bins <- function(df) {
        new_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
        colnames(new_df) <- colnames(df)
        
        i <- 1
        while (i <= nrow(df)) {
            current_bin_size <- df[[bin_size_col]][i]
            
            # If the bin size of the first bin is too small, combine with next bin
            if (i == 1 && current_bin_size < mean_size - sd_size) {
                if (verbose) {
                    message(paste("Combining first bin with bin", i+1))
                }
                df[i+1, l_bound_col] <- df[i, l_bound_col]
                df[i+1, bin_size_col] <- df[i, bin_size_col] + df[i+1, bin_size_col]
                i <- i + 2
            }
            # If bin size is too small (not the first bin), combine with previous bin
            else if (current_bin_size < mean_size - sd_size) {
                if (verbose) {
                    message(paste("Combining bin", i, "with previous bin"))
                }
                new_df[nrow(new_df), u_bound_col] <- df[i, u_bound_col]
                new_df[nrow(new_df), bin_size_col] <- new_df[nrow(new_df), bin_size_col] + current_bin_size
                i <- i + 1
            } 
            # If bin size is too large, split it
            else if (current_bin_size > mean_size + sd_size) {
                num_splits <- round(current_bin_size / mean_size)
                split_sizes <- rep(floor(current_bin_size / num_splits), num_splits)
                remaining <- current_bin_size - sum(split_sizes)
                split_sizes[1:remaining] <- split_sizes[1:remaining] + 1
                
                for (k in 1:num_splits) {
                    l_bound <- df[i, l_bound_col] + (sum(split_sizes[1:k-1]) / current_bin_size) * (df[i, u_bound_col] - df[i, l_bound_col])
                    u_bound <- l_bound + (split_sizes[k] / current_bin_size) * (df[i, u_bound_col] - df[i, l_bound_col])
                    new_row <- data.frame(
                        l_bound_col = l_bound,
                        u_bound_col = u_bound,
                        bin_size_col = split_sizes[k],
                        binPTime_col = NA,
                        stringsAsFactors = FALSE
                    )
                    colnames(new_row) <- c(l_bound_col, u_bound_col, bin_size_col, binPTime_col)
                    new_df <- rbind(new_df, new_row)
                }
                if (verbose) {
                    message(paste("Splitting bin", i, "into", num_splits, "bins"))
                }
                i <- i + 1
            } 
            # Otherwise, just keep the bin
            else {
                new_df <- rbind(new_df, df[i, ])
                i <- i + 1
            }
        }
        
        # Reset the pseudotime column
        new_df[[binPTime_col]] <- 1:nrow(new_df)
        
        return(new_df)
    }
    
    # Optimize bins iteratively until there's no change
    iteration_count <- 0
    previous_df <- df
    while (TRUE) {
        iteration_count <- iteration_count + 1
        if (verbose) {
            message(paste("Optimization iteration:", iteration_count))
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
    
    return(optimized_df)
}
