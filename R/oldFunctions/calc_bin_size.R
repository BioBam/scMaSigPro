calc_bin_size <- function(x) {
  # Break and calc length
  size <- length(
    c(stringr::str_split(x[["cluster.members"]], "\\|"))[[1]]
  )

  return(as.numeric(size))
}
