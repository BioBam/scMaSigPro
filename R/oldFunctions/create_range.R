create_range <- function(x) {
  # Convert the factor into Character
  y <- as.character(x[["bin"]])
  y <- y %>% stringr::str_remove_all(pattern = "\\[|\\]|\\(|\\)")
  y1 <- as.numeric(sapply(strsplit(y, ","), "[", 1))
  y2 <- as.numeric(sapply(strsplit(y, ","), "[", 2))
  rangeVec <- c(y1, y2, x[["bin_size"]], x[["customTime"]])
  return(as.numeric(rangeVec))
}
