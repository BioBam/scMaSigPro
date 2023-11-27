getPattern <- function(frame, trend, group, groups_vector) {
  # Get path related columns
  frame <- frame[, groups_vector == group, drop = FALSE]

  # Drop beta0 if present
  if ("beta0" %in% colnames(frame)) {
    frame <- frame[, colnames(frame) != "beta0", drop = FALSE]
  }

  # Subset by trend
  if (trend == "up") {
    frame <- frame[rowSums(frame) > 0, , drop = FALSE]
  } else if (trend == "down") {
    frame <- frame[rowSums(frame) <= 0, , drop = FALSE]
  }

  return(frame)
}
