#' @title Polynomial term selection methods from maSigPro
#'
#' @description
#' All function are modified to include the offset
#'
#' @references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
#' Bioinformatics 22, 1096-1102}
#'
#' @author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@@ua.es}}
#'
#' @keywords internal
# Two Way Step-Back
sc.two.ways.stepback <- function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001, useOffset, useWeight, max_it) {
  OUT <- NULL
  lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
  result <- summary(lm1)$coefficients[, 4]
  max <- max(result[-1], na.rm = TRUE)
  d <- d[, names(result)[-1]]
  while (max > alfa) {
    varout <- names(result)[result == max]
    # Clean String
    varout <- clean_string(varout, action = "remove")

    pos <- maSigPro::position(matrix = d, vari = varout)
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
      lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
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
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
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
    lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
    result <- summary(lm1)$coefficients[, 4]
    max <- max(result[-1], na.rm = TRUE)
    if (length(result[-1]) == 1) {
      max <- result[-1]
      if (max > alfa) {
        max <- 0
        lm1 <- glm(y ~ 1, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
      }
    }
  }
  return(lm1)
}

# Two-way step forward
sc.two.ways.stepfor <-
  function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001, useOffset, useWeight, weights = useWeight, max_it) {
    pval <- NULL
    design <- NULL
    j <- 1
    resul0 <- summary(glm(y ~ ., data = d, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it))$coefficients[, 4]
    d <- as.data.frame(d[, names(resul0)[-1]])
    for (i in 1:ncol(d)) {
      sub <- cbind(design, d[, i])
      sub <- as.data.frame(sub)
      lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
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
      result2 <- summary(glm(y ~ ., data = design, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it))$coefficients[
        ,
        4
      ]
      max <- max(result2[-1], na.rm = TRUE)
      while (max > alfa) {
        varout <- names(result2)[result2 == max]
        # Clean String
        varout <- clean_string(varout, action = "remove")
        pos <- maSigPro::position(matrix = design, vari = varout)
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
        result2 <- summary(glm(y ~ ., data = design, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it))$coefficients[
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
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
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
      lm1 <- glm(y ~ 1, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
    } else {
      lm1 <- glm(y ~ ., data = design, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
    }
    return(lm1)
  }

# Stepback
sc.stepback <- function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001, useOffset, useWeight, max_it) {
  lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
  result <- summary(lm1)
  max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
  if (length(result$coefficients[, 4][-1]) == 1) {
    if (max > alfa) {
      max <- 0
      lm1 <- glm(y ~ 1, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
    }
  }
  while (max > alfa) {
    varout <- as.character(names(result$coefficients[, 4])[result$coefficients[
      ,
      4
    ] == max][1])

    # Clean String
    varout <- clean_string(varout, action = "remove")

    pos <- maSigPro::position(matrix = d, vari = paste(varout))
    d <- d[, -pos]
    if (length(result$coefficients[, 4][-1]) == 2) {
      min <- min(result$coefficients[, 4][-1], na.rm = TRUE)
      lastname <- names(result$coefficients[, 4][-1])[result$coefficients[, 4][-1] == min]
    }
    if (is.null(dim(d))) {
      d <- as.data.frame(d)
      colnames(d) <- lastname
    }
    lm1 <- glm(y ~ ., data = d, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
    result <- summary(lm1)
    max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
    if (length(result$coefficients[, 4][-1]) == 1) {
      max <- result$coefficients[, 4][-1]
      if (max > alfa) {
        max <- 0
        lm1 <- glm(y ~ 1, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
      }
    }
  }
  return(lm1)
}

# Stepfor
sc.stepfor <- function(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 0.00001, useOffset, useWeight, max_it) {
  pval <- NULL
  design <- NULL
  j <- 1
  resul0 <- summary(glm(y ~ ., data = d, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it))$coefficients[, 4]
  d <- as.data.frame(d[, names(resul0)[-1]])
  for (i in 1:ncol(d)) {
    sub <- cbind(design, d[, i])
    sub <- as.data.frame(sub)
    lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
    result <- summary(lm2)
    pval[i] <- result$coefficients[, 4][j + 1]
  }
  min <- min(pval)
  while (min < alfa) {
    b <- pval == min
    c <- c(1:length(pval))
    pos <- c[b]
    pos <- pos[!is.na(pos)][1]
    design <- cbind(design, d[, pos])
    design <- as.data.frame(design)
    colnames(design)[j] <- colnames(d)[pos]
    j <- j + 1

    if (ncol(d) == 2) {
      lastname <- colnames(d)[!b]
    }
    d <- as.data.frame(d[, -pos])
    if (ncol(d) == 1) {
      colnames(d) <- lastname
    }

    pval <- NULL
    if (ncol(d) != 0) {
      for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub)
        lm2 <- glm(y ~ ., data = sub, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1]
      }
      min <- min(pval, na.rm = TRUE)
    } else {
      min <- 1
    }
  }
  if (is.null(design)) {
    lm1 <- glm(y ~ 1, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
  } else {
    lm1 <- glm(y ~ ., data = design, family = family, epsilon = epsilon, offset = useOffset, weights = useWeight, maxit = max_it)
  }
  return(lm1)
}
