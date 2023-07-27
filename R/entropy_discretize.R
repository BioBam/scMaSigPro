#' entropy_discretize
#'
#' @description
#' This function discretizes a continuous time series column into bins of equal size using entropy-based binning method. It automatically calculates the optimal number of bins using one of the supported methods.
#' The bin sizes are also calculated and merged with the input design_table.
#'
#' @param design_table A data.frame. The dataset containing the time series column to be discretized.
#' @param time_col A character string. The name of the column in design_table that contains the time series data to be discretized.
#' @param method A character string (default = "Sturges"). The method to be used to estimate the optimal number of bins. Currently, this function supports "Sturges" method.
#' @param drop.fac A numeric value (default = 0.5). The factor by which to decrease the number of bins if the initial binning results in too many empty bins. The optimal number of bins is recalculated until this criteria is met.
#' @param verbose A boolean (default = TRUE). If TRUE, detailed messages about the process (e.g. number of bins calculated) will be printed.
#'
#' @return
#' A data.frame that contains the original data plus additional columns:
#' - 'bin' : the bin number
#' - 'bin_size' : the size of the bin
#' - 'binned_time' : the interval range of each bin
#' This function returns the merged data.frame with new discretized time_col, preserving the original rownames.
#'
#' @details
#' This function performs the following steps:
#' - Adds a new column 'cell' to the input data.frame which copies the row names.
#' - Extracts the time series data from the specified column of the input data.frame.
#' - Calculates the optimal number of bins using the specified method.
#' - Prints the estimated number of bins if verbose is set to TRUE.
#' - Discretizes the time series data into bins using the entropy-based binning method.
#' - Merges the original data.frame with the new binned time series data.
#' - Removes the 'cell' column and sets the row names back to the original row names of the input data.frame.
#' - Returns the merged data.frame.
#'
#' @examples
#' \dontrun{
#' entropy_discretize(design_table = data.frame, time_col = "time",
#' method = "Sturges", drop.fac = 0.5, verbose = TRUE)
#' }
#' 
#' @seealso \code{\link{estBinSize}}, \code{\link{discretize}}, \code{\link{create_range}}
#' 
#' @export

entropy_discretize <- function(design_table, time_col,
                               method = "Sturges",
                               drop.fac = 0.5,
                               verbose = TRUE) {

  # Add a column
  design_table$cell <- rownames(design_table)

  # Extract the time information as a vector
  time_vector <- design_table[, time_col]
  length_n <- length(time_vector)
  
  # Calculate Optimal Number of Bins
  estBins <- estBinSize(time_vector = time_vector, nPoints = length_n,
                        drop_fac = drop.fac, method = method)
 
  # Client-Verbose
  if (verbose){message(paste("Estimated Bin Sizes =", estBins, "with",
                    method, "binning for", length_n, "time points."))}
  
  # Calculate Bin intervals with entropy
  bin_intervals <- as.data.frame(discretize(time_vector, numBins = estBins, r = range(time_vector)))
  
  # Clean the table before merge
  colnames(bin_intervals) <- c("bin", "bin_size")
  bin_intervals$binned_time <- rownames(bin_intervals)

  # Create the bin table
  bin_table <- as.data.frame(t(as.data.frame(apply(bin_intervals, 1, create_range))))
  colnames(bin_table) <- c("from", "to", "bin_size", "binnedTime")

  # Merge with design table
  processed_design_table <- as.data.frame(left_join(design_table, bin_table,
    by = join_by(closest(!!time_col >= from), closest(!!time_col <= to))
  ))

  # Remove cell column
  rownames(processed_design_table) <- processed_design_table$cell

  # Drop cell column
  processed_design_table <- processed_design_table[, colnames(processed_design_table) != "cell"]
  
  return(processed_design_table)
}
