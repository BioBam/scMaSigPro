"two.ways.stepfor" <-
  function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001) {
    pval <- NULL
    design <- NULL
    j <- 1
    resul0 <- summary(glm(y ~ ., data = d, family = family, epsilon = epsilon))$coefficients[, 4]
    d <- as.data.frame(d[, names(resul0)[-1]])
    for (i in 1:ncol(d)) {
      sub <- cbind(design, d[, i])
      sub <- as.data.frame(sub)
      lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon)
      result <- summary(lm2)
      pval[i] <- result$coefficients[, 4][j + 1]
    }
    min <- min(pval, na.rm = TRUE)
    while (min <= alfa) {
      b <- pval == min
      c <- c(1:length(pval))
      pos <- c[b]
      pos <- pos[!is.na(pos)][1]
      design <- cbind(design, d[, pos])
      design <- as.data.frame(design)
      colnames(design)[j] <- colnames(d)[pos]
      if (ncol(d) == 2) {
        lastname <- colnames(d)[colnames(d) != colnames(d)[pos]]
      }
      d <- d[, -pos]
      if (is.null(dim(d))) {
        d <- as.data.frame(d)
        colnames(d) <- lastname
      }
      result2 <- summary(glm(y ~ ., data = design, family = family, epsilon = epsilon))$coefficients[
        ,
        4
      ]
      max <- max(result2[-1], na.rm = TRUE)
      while (max > alfa) {
        varout <- names(result2)[result2 == max]
        pos <- position(matrix = design, vari = varout)
        d <- as.data.frame(cbind(d, design[, pos]))
        x <- ncol(d)
        colnames(d)[x] <- colnames(design)[pos]
        if (ncol(design) == 2) {
          min <- min(result2[-1], na.rm = TRUE)
          lastname <- names(result2)[result2 == min]
        }
        design <- design[, -pos]
        if (is.null(dim(design))) {
          design <- as.data.frame(design)
          colnames(design) <- lastname
        }
        result2 <- summary(glm(y ~ ., data = design, family = family, epsilon = epsilon))$coefficients[
          ,
          4
        ]
        max <- max(result2[-1], na.rm = TRUE)
      }
      j <- ncol(design) + 1
      pval <- NULL
      for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub)
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1]
      }
      min <- min(pval, na.rm = TRUE)
      if (ncol(d) == 1) {
        if (min <= alfa) {
          design <- cbind(design, d[, 1])
          design <- as.data.frame(design)
          colnames(design)[j] <- colnames(d)[1]
        }
        min <- 1
      }
    }
    if (is.null(design)) {
      lm1 <- glm(y ~ 1, family = family, epsilon = epsilon)
    } else {
      lm1 <- glm(y ~ ., data = design, family = family, epsilon = epsilon)
    }
    return(lm1)
  }
