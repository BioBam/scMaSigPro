entropy_discretize <- function(design_table, time_col,
                               method = "Sturges", drop.fac = 0.5) {
  require(tidyverse)

  # Add a column
  design_table$cell <- rownames(design_table)

  # Extract the time information as a vector
  time_vector <- design_table[, time_col]
  length_n <- length(time_vector)

  # Calculate optimum Number of bins
  if (method == "Freedman.Diaconis") {
    estBins <- 2 * IQR(time_vector) / length_n^(1 / 3)
    estBins <- drop.fac * estBins
  } else if (method == "Sqrt") {
    estBins <- length_n^(1 / 2)
    estBins <- drop.fac * estBins
  } else if (method == "Sturges") {
    estBins <- log2(length_n) + 1
    estBins <- drop.fac * estBins
  } else if (method == "Rice") {
    estBins <- 2 * length_n^(1 / 3)
    estBins <- drop.fac * estBins
  } else if (method == "Doane") {
    sigma <- ((6 * (length_n - 2)) / ((length_n + 1) * (length_n + 3)))^(1 / 2)
    sk <- moments::skewness(time_vector)
    estBins <- 1 + log2(length_n) + log2(1 + (abs(sk) / sigma))
    estBins <- drop.fac * estBins
  } else if (method == "Scott.Normal") {
    estBins <- 3.49 * abs(sd(time_vector)) / length_n^(1 / 3)
    estBins <- drop.fac * estBins
  }

  # Calculate Bin intervals with entropy
  suppressPackageStartupMessages(require(entropy))
  bin_intervals <- as.data.frame(discretize(time_vector, numBins = estBins, r = range(time_vector)))

  # Clean the table before merge
  colnames(bin_intervals) <- c("bin", "bin_size")
  bin_intervals$customTime <- rownames(bin_intervals)

  # Create the bin table
  bin_table <- as.data.frame(t(as.data.frame(apply(bin_intervals, 1, create_range))))
  colnames(bin_table) <- c("from", "to", "bin_size", "binnedTime")

  # Merge with design table
  processed_design_table <- as.data.frame(left_join(design_table, bin_table,
    by = join_by(closest(!!time_col >= from), closest(!!time_col <= to))
  )) # %>%
  # select(names(design_table), binnedTime))


  # Remove cell column
  rownames(processed_design_table) <- processed_design_table$cell

  # Drop cell column
  processed_design_table <- processed_design_table[, colnames(processed_design_table) != "cell"]

  return(processed_design_table)
}
