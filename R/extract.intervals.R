#' @title Extract Intervals
#' 
#' @param time.vector 
#' @param nBins 
#' @param bin 
#' @param bin.size 
#' @param lbond
#' @param ubond 
#' 
#' @keywords internal

extract.intervals <- function(time.vector, nBins = 1, bin, bin.size, lbound, ubound) {
  # Create Dataframe
  new_range_current <- as.data.frame(entropy::discretize(
    time.vector,
    numBins = nBins, r = range(time.vector)
  ))
  # Set columns
  colnames(new_range_current) <- c(bin, bin.size)

  # Current bin new table
  new_bin_table_current <- as.data.frame(
    t(as.data.frame(
      apply(
        new_range_current, 1, create_range,
        bin_size_colname = bin.size,
        bin_colname = bin,
        verbose = FALSE
      )
    ))
  )
  # Set column names
  colnames(new_bin_table_current) <- c(lbound, ubound, bin.size)

  # Return
  return(new_bin_table_current)
}
