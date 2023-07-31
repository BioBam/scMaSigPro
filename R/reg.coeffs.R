"reg.coeffs" <-
  function(
      coefficients, indepen = groups.vector[nchar(groups.vector) == min(nchar(groups.vector))][1], groups.vector,
      group) {
    coefficients <- as.numeric(coefficients)
    c0 <- coefficients[1]
    if (length(coefficients) == length(groups.vector)) groups.vector <- groups.vector[-1]
    coefficients <- coefficients[-1]
    c1 <- coefficients[groups.vector == indepen]
    c01 <- c(c0, c1)
    if (group != indepen) {
      group <- paste(group, "vs", indepen, sep = "")
      c2 <- coefficients[group == groups.vector]
      if (length(c2) < length(c01)) {
        c2 <- c(0, c2)
      } else if (length(c2) != length(c01)) {
        stop("incomplete coeffcients")
      }
      reg.coeff <- c01 + c2
    } else {
      reg.coeff <- c01
    }
    reg.coeff
  }
