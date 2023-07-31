IsoModel <- function(
    data, gen, design = NULL, Q = 0.05, min.obs = 6, minorFoldfilter = NULL,
    counts = FALSE, family = NULL, theta = 10, epsilon = 1e-05) {
  #---------------------------------------------------------------------------------------------------
  # data is a matrix containing isoform expression. Isoforms must be in rows and experimental conditions in columns
  # gen is a vector with the name of the gene each isoform belongs to
  #---------------------------------------------------------------------------------------------------

  Genes <- unique(gen)
  g <- length(Genes)

  if (is.null(family)) {
    if (!counts) {
      family <- gaussian()
    }
    if (counts) {
      family <- negative.binomial(theta)
    }
  }

  #---------------------------------------------------------------------------------------------------
  # STEP -1: Remove cases with low expressed isoforms:
  #---------------------------------------------------------------------------------------------------

  print(paste(nrow(data), "transcripts"))
  print(paste(length(unique(gen)), "genes"))

  if (!is.null(minorFoldfilter)) {
    print("Removing low expressed minor isoforms")
    moreOne <- names(which(table(gen) > 1))
    iso.sel <- NULL
    gene.sel <- NULL
    for (i in moreOne) {
      which(gen == i)
      gene.data <- data[which(gen == i), ]
      isoSUM <- apply(gene.data, 1, sum)
      major <- names(which(isoSUM == max(isoSUM)))[1]
      minors <- names(which(isoSUM != max(isoSUM)))
      div <- as.numeric(matrix(rep(gene.data[major, ], length(minors)), ncol = ncol(data), length(minors), byrow = T)) / as.matrix(gene.data[minors, ])
      is <- names(which(apply(div, 1, min, na.rm = T) < minorFoldfilter))
      iso.sel <- c(iso.sel, is)
      gene.sel <- c(gene.sel, rep(i, length(is)))
    }
    data <- data[iso.sel, ]
    gen <- gene.sel
    print(dim(data))
    print(length(gen))
    print("Done")
    print(paste(nrow(data), "remaining transcripts"))
    print(paste(length(unique(gen)), "remaining genes"))
  }

  #---------------------------------------------------------------------------------------------------
  #  STEP 0: Remove cases with 1 transcript:
  #---------------------------------------------------------------------------------------------------

  NT <- tapply(rownames(data), gen, length)
  Genes1 <- names(which(NT != 1))
  data1 <- data[gen %in% Genes1, ]
  gen1 <- gen[gen %in% Genes1]
  Genes1 <- unique(gen1)
  g <- length(Genes1)
  dis <- as.data.frame(design$dis)
  mycolnames <- colnames(dis)

  #---------------------------------------------------------------------------------------------------
  # STEP 1: Gene models comparison.
  #---------------------------------------------------------------------------------------------------

  pval <- NULL

  for (i in 1:g)
  {
    div <- c(1:round(g / 100)) * 100
    if (is.element(i, div)) {
      print(paste(c("fitting gene", i, "out of", g), collapse = " "))
    }

    zz <- data1[gen1 == Genes1[i], ]
    nt <- nrow(zz)

    dis.gen <- REP(dis, nt)
    y <- c(t(as.matrix(zz)))
    transcript <- factor(rep(c(1:nt), each = ncol(zz)))
    ydis <- cbind(y, dis.gen, transcript)

    model0 <- glm(Formula0(mycolnames), data = ydis, family = family, epsilon = epsilon)
    model1 <- glm(Formula1(mycolnames), data = ydis, family = family, epsilon = epsilon)

    if (family$family == "gaussian") {
      pvali <- anova(model0, model1, test = "F")[2, 6]
    } else {
      pvali <- anova(model0, model1, test = "Chisq")[2, 5]
    }
    names(pvali) <- Genes1[i]
    pval <- c(pval, pvali)
  }
  num.genes <- sum(p.adjust(pval) < Q, na.rm = TRUE)
  selected.genes <- names(sort(p.adjust(pval))[1:num.genes])

  #---------------------------------------------------------------------------------------------------
  # STEP 2: p.vector and T.fit to the transcripts that belong to selected.genes
  #---------------------------------------------------------------------------------------------------
  data2 <- data[gen %in% selected.genes, ]
  gen2 <- gen[gen %in% selected.genes]
  pvector2 <- p.vector(data2, design, counts = counts, item = "isoform")
  Tfit2 <- T.fit(pvector2, item = "isoform")

  #---------------------------------------------------------------------------------------------------
  # Output
  #---------------------------------------------------------------------------------------------------

  ISO.SOL <- list(data, gen, design, selected.genes, pvector2, Tfit2)
  names(ISO.SOL) <- c("data", "gen", "design", "DSG", "pvector", "Tfit")
  ISO.SOL
}


#---------------------------------------------------------------------------------------------------
# Auxiliar internal functions: REP, Formula0, Formula1
#---------------------------------------------------------------------------------------------------

REP <- function(D, k) {
  r <- nrow(D)
  c <- ncol(D)
  DD <- NULL
  for (i in 1:c)
  {
    DDi <- rep(D[, i], k)
    DD <- cbind(DD, DDi)
  }
  colnames(DD) <- colnames(D)
  as.data.frame(DD)
}

#---------------------------------------------------------------------------

Formula0 <- function(names) {
  formula <- "y~"

  if (length(names) == 1) {
    formula <- paste(formula, names[1], "+ transcript")
  } else if (length(names) > 1) {
    for (i in 1:(length(names)))
    {
      formula <- paste(formula, names[i], "+")
    }
    formula <- paste(formula, "transcript")
  }
  formula <- as.formula(formula)
  formula
}

#---------------------------------------------------------------------------

Formula1 <- function(names) {
  formula <- "y~"

  if (length(names) == 1) {
    formula <- paste(formula, names[1], "* transcript")
  } else if (length(names) > 1) {
    formula <- paste(formula, "(")
    for (i in 1:(length(names) - 1))
    {
      formula <- paste(formula, names[i], "+")
    }
    formula <- paste(formula, names[length(names)])
    formula <- paste(formula, ") * transcript")
  }
  formula <- as.formula(formula)
  formula
}
