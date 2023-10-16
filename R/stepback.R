stepback <- function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001) {
  lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon)
  result <- summary(lm1)
  max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
  if (length(result$coefficients[, 4][-1]) == 1) {
    if (max > alfa) {
      max <- 0
      lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
    }
  }
  while (max > alfa) {
    varout <- names(result$coefficients[, 4])[result$coefficients[
      ,
      4
    ] == max][1]
    pos <- position(matrix = d, vari = varout)
    d <- d[, -pos]
    if (length(result$coefficients[, 4][-1]) == 2) {
      min <- min(result$coefficients[, 4][-1], na.rm = TRUE)
      lastname <- names(result$coefficients[, 4][-1])[result$coefficients[, 4][-1] == min]
    }
    if (is.null(dim(d))) {
      d <- as.data.frame(d)
      colnames(d) <- lastname
    }
    lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon)
    result <- summary(lm1)
    max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
    if (length(result$coefficients[, 4][-1]) == 1) {
      max <- result$coefficients[, 4][-1]
      if (max > alfa) {
        max <- 0
        lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
      }
    }
  }
  return(lm1)
}
