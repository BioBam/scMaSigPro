extract.group.components <- function(coeff.vec, groups) {
  # Convert NA to zero
  coeff.vec[is.na(coeff.vec)] <- 0

  # Extract A
  A <- coeff.vec[c(T, F)]

  print(A)

  # Extract B
  B <- A + coeff.vec[c(F, T)]

  print(B)

  stop()

  # Names
  names(A) <- paste(groups[1], "beta", c(0:(length(A) - 1)), sep = "_")
  names(B) <- paste(groups[2], "beta", c(0:(length(B) - 1)), sep = "_")

  # Cbind
  group.coeff <- c(A, B)

  # Return
  return(group.coeff)
}
