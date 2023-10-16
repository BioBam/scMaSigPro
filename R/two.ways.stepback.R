"two.ways.stepback" <-
  function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001) {
    OUT <- NULL
    lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon)
    result <- summary(lm1)$coefficients[, 4]
    max <- max(result[-1], na.rm = TRUE)
    d <- d[, names(result)[-1]]
    while (max > alfa) {
      varout <- names(result)[result == max]
      pos <- position(matrix = d, vari = varout)
      OUT <- as.data.frame(cbind(OUT, d[, pos]))
      x <- ncol(OUT)
      colnames(OUT)[x] <- colnames(d)[pos]
      if (ncol(d) == 2) {
        min <- min(result[-1], na.rm = TRUE)
        lastname <- names(result)[result == min]
      }
      d <- d[, -pos]
      if (is.null(dim(d))) {
        d <- as.data.frame(d)
        colnames(d) <- lastname
      }
      j <- ncol(d) + 1
      pval <- NULL
      for (i in 1:ncol(OUT)) {
        sub <- cbind(d, OUT[, i])
        sub <- as.data.frame(sub)
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon)
        result <- summary(lm2)$coefficients[, 4]
        pval[i] <- result[j + 1]
      }
      min <- min(pval, na.rm = TRUE)
      while (min <= alfa) {
        b <- pval == min
        c <- c(1:length(pval))
        pos <- c[b]
        d <- cbind(d, OUT[, pos])
        d <- as.data.frame(d)
        colnames(d)[j] <- colnames(OUT)[pos]
        if (ncol(OUT) == 2) {
          max <- max(pval, na.rm = TRUE)
          b <- pval == max
          c <- c(1:length(pval))
          last <- c[b]
          lastname <- colnames(OUT)[last]
        }
        OUT <- OUT[, -pos]
        if (is.null(dim(OUT))) {
          OUT <- as.data.frame(OUT)
          colnames(OUT) <- lastname
        }
        j <- ncol(d) + 1
        pval <- NULL
        for (i in 1:ncol(OUT)) {
          sub <- cbind(d, OUT[, i])
          sub <- as.data.frame(sub)
          lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon)
          result <- summary(lm2)
          pval[i] <- result$coefficients[, 4][j + 1]
        }
        min <- min(pval, na.rm = TRUE)
        if (ncol(OUT) == 1) {
          if (min <= alfa) {
            d <- cbind(d, OUT[, 1])
            d <- as.data.frame(d)
            colnames(d)[j] <- colnames(OUT)[1]
          }
          min <- 1
        }
      }
      lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon)
      result <- summary(lm1)$coefficients[, 4]
      max <- max(result[-1], na.rm = TRUE)
      if (length(result[-1]) == 1) {
        max <- result[-1]
        if (max > alfa) {
          max <- 0
          lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
        }
      }
    }
    return(lm1)
  }
