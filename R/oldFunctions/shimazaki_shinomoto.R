shin_shin <- function(design_table, time_col = "Step",
                      min_bin = 16, max_bin = 100) {
  # Add a column
  design_table$cell <- rownames(design_table)

  # Extract the time information as a vector
  time_vector <- design_table[, time_col]
  # time_vector <- somevector

  # Estimate the minimum and max
  min_time <- min(time_vector)
  max_time <- max(time_vector)

  # Generate a vector of possible bin-widths
  bin.vec <- seq(min_bin, max_bin - 1, 1)

  # Bin-size vector
  bin.size.vec <- (max_time - min_time) / bin.vec

  # zero vector
  zero.vec <- c()

  for (i in c(1:length(bin.vec))) {
    edge <- seq(from = min_time, to = max_time, length.out = bin.vec[i] + 1)
    h <- hist(time_vector, breaks = edge)
    h.counts <- h$counts
    h.counts.mean <- mean(h.counts)
    h.v <- sum((h.counts - h.counts.mean)**2) / bin.vec[i]
    val <- (2 * h.counts.mean - h.v) / ((bin.size.vec[i])**2)

    zero.vec <- c(zero.vec, val)
  }

  zero.vec <- zero.vec[!is.na(zero.vec)]


  min_val <- min(zero.vec)

  print(which(min_val == zero.vec))

  print(zero.vec)
  # optD <- bin.size.vec[which(min(zero.vec) == zero.vec)]
  optBin <- bin.vec[which(min_val == zero.vec)]
  edge <- seq(
    from = min_time, to = max_time,
    length.out = optBin
  )

  print(edge)
}
