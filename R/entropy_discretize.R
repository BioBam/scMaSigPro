#' entropy_discretize
#'
#' @description
#' This function discretizes a continuous time series column into bins of equal size using entropy-based binning method. It automatically calculates the optimal number of bins using one of the supported methods.
#' The bin sizes are also calculated and merged with the input design_table.
#'
#' @param design_table A data.frame. The dataset containing the time series column to be discretized.
#' @param time_col A character string. The name of the column in design_table that contains the time series data to be discretized.
#' @param path_col A character string. The name of the column in design_table that contains the Lineage information.
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
#' entropy_discretize(
#'   design_table = data.frame, time_col = "time",
#'   method = "Sturges", drop.fac = 0.5, verbose = TRUE
#' )
#' }
#'
#' @seealso \code{\link{estBinSize}}, \code{\link{discretize}}, \code{\link{create_range}}
#'
#' @export

entropy_discretize <- function(design_table, time_col,
                               path_col,
                               method = "Sturges",
                               drop.fac = 0.5,
                               verbose = TRUE) {
  # Checks
  assert_that(time_col %in% colnames(design_table),
    msg = paste0("'", time_col, "' does not exist in design_table")
  )
  assert_that(path_col %in% colnames(design_table),
    msg = paste0("'", path_col, "' does not exist in design_table")
  )
  assert_that(drop.fac >= 0.3 & drop.fac <= 1,
    msg = "Invalid value for 'drop.fac'. It should be between 0.3 and 1."
  )

  # Add a column
  design_table$cell <- rownames(design_table)

  # Get the avaible paths
  avail.paths <- as.vector(unique(design_table[[path_col]]))

  # Check for path
  assert_that(length(avail.paths) >= 2,
    msg = "Invalid number of paths detected. Please make sure that dataset has atleast two paths"
  )
  
  # Determine the number of cores
  num_cores <- detectCores() - 1
  
  # Apply transformations on data
  #discrete.list <- mclapply(avail.paths, function(path, design.frame = design_table,
    discrete.list <- lapply(avail.paths, function(path, design.frame = design_table,
                                                  drop_fac = drop.fac, path.col = path_col,
                                                  time.col = time_col, method.bin = method) {
      # Get the cells belonging to path
      path.frame <- design.frame[design.frame[[path.col]] == path, , drop = F]
      
      # Extract the time information as a vector
      time_vector <- path.frame[, time.col]
      length_n <- length(time_vector)
      
      # Validation
      if (length_n <= 7) {
          message(paste("Time points are already less than 7 in", path))
      }
      
      # Calculate Optimal Number of Bins
      tryCatch(
          expr = {
              estBins <- estBinSize(
                  time_vector = time_vector, nPoints = length_n,
                  drop_fac = drop.fac, method = method.bin
              )
              
              if (verbose) {
                  message(paste(
                      "Estimated Bin Sizes =", estBins, "with",
                      method, "binning for", length_n, "time points for", path
                  ))
              }
          },
          error = function(e) {
              message(paste("Error message: ", e$message))
              stop("Unable to estimate bin size")
          }
      )
      
      # Calculate Bin intervals with entropy
      bin_intervals <- as.data.frame(discretize(time_vector, numBins = estBins, r = range(time_vector)))
      
      # Client-Verbose
      if (verbose) {
          message(paste(
              "For", path, ",", length_n, "time points has been compressed to", nrow(bin_intervals), "bins"
          ))
      }
      
      # Clean the table before merge
      colnames(bin_intervals) <- c("bin", "bin_size")
      bin_intervals$binned_time <- rownames(bin_intervals)
      
      # Create the bin table
      bin_table <- as.data.frame(t(as.data.frame(apply(bin_intervals, 1, create_range))))
      colnames(bin_table) <- c("from", "to", "bin_size", "binnedTime")
      
      # Merge with design table
      processed_design_table <- as.data.frame(left_join(path.frame, bin_table,
                                                        by = join_by(closest(!!time.col >= from), closest(!!time.col <= to))
      ))
      
      return(processed_design_table)
  #}, mc.cores = num_cores)
  })
  
  # Bind rows
  discrete.frame <- bind_rows(discrete.list)
  
  # Remove cell column and set rows
  processed_design_table <- discrete.frame %>%
      rownames_to_column(var = "row_id") %>%
      mutate(row_id = design_table$cell) %>%
      column_to_rownames(var = "row_id") %>% as.data.frame()
  
  # Drop cell
  processed_design_table <- processed_design_table[, colnames(processed_design_table) != "cell"]
  
  # Convert the resulting tibble to a dataframe
  processed_design_table <- as.data.frame(processed_design_table)
  
  return(processed_design_table)
  
}
