p.vector <- function(data, design, Q = 0.05, MT.adjust = "BH", min.obs = 6, counts = FALSE, family = NULL, theta = 10, epsilon = 0.00001, item = "gene") {
  if (is.data.frame(design) || is.matrix(design)) {
    dis <- design
    groups.vector <- NULL
    edesign <- NULL
  } else if (is.list(design)) {
    dis <- as.data.frame(design$dis)
    groups.vector <- design$groups.vector
    edesign <- design$edesign
  }

  if (is.null(family)) {
    if (!counts) {
      family <- gaussian()
    }
    if (counts) {
      family <- negative.binomial(theta)
    }
  }

  dat <- as.matrix(data)
  dat <- dat[, as.character(rownames(dis))]
  G <- nrow(dat)

  # Removing rows with many missings:

  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min.obs, ]

  # Removing rows with all ceros:

  sumatot <- apply(dat, 1, sum)
  counts0 <- which(sumatot == 0)
  if (length(counts0) > 0) {
    dat <- dat[-counts0, ]
  }

  g <- dim(dat)[1]
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  p.vector <- vector(mode = "numeric", length = g)

  #----------------------------------------------------------------------
  for (i in 1:g) {
    y <- as.numeric(dat[i, ])

    div <- c(1:round(g / 100)) * 100
    if (is.element(i, div)) {
      print(paste(c("fitting ", item, i, "out of", g), collapse = " "))
    }

    model.glm <- glm(y ~ ., data = dis, family = family, epsilon = epsilon)
    if (model.glm$null.deviance == 0) {
      p.vector[i] <- 1
    } else {
      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)

      if (family$family == "gaussian") {
        test <- anova(model.glm.0, model.glm, test = "F")
        if (is.na(test[6][2, 1])) {
          p.vector[i] <- 1
        } else {
          p.vector[i] <- test[6][2, 1]
        }
      } else {
        test <- anova(model.glm.0, model.glm, test = "Chisq")
        if (is.na(test[5][2, 1])) {
          p.vector[i] <- 1
        } else {
          p.vector[i] <- test[5][2, 1]
        }
      }
    }
  }

  #----------------------------------------------------------------------
  p.adjusted <- p.adjust(p.vector, method = MT.adjust, n = length(p.vector))
  genes.selected <- rownames(dat)[which(p.adjusted <= Q)]

  FDR <- sort(p.vector)[length(genes.selected)]

  SELEC <- as.matrix(as.data.frame(dat)[genes.selected, ])
  if (nrow(SELEC) == 0) {
    print("no significant genes")
  }

  p.vector <- as.matrix(p.vector)
  rownames(p.vector) <- rownames(dat)
  colnames(p.vector) <- c("p.value")

  #-------------------------------------------------------------------------

  output <- list(
    SELEC, p.vector, p.adjusted, G, g, FDR,
    nrow(SELEC), dis, dat, min.obs, Q, groups.vector, edesign, family
  )
  names(output) <- c(
    "SELEC", "p.vector", "p.adjusted", "G",
    "g", "FDR", "i", "dis", "dat", "min.obs", "Q", "groups.vector",
    "edesign", "family"
  )
  output
}
