#' Make regression fit for time series gene expression experiments
#'
#' \code{sc.p.vector} performs a regression fit for each gene taking all variables present in the model given by a regression matrix
#' and returns a list of FDR corrected significant genes.
#'
#' @param scmpObj matrix containing normalized gene expression data. Genes must be in rows and arrays in columns.
#' @param Q significance level. Default is 0.05.
#' @param MT.adjust argument to pass to \code{p.adjust} function indicating the method for multiple testing adjustment of p.value.
#' @param min.obs genes with less than this number of true numerical values will be excluded from the analysis.
#'   Minimum value to estimate the model is (degree+1) x Groups + 1. Default is 6.
#' @param counts a logical indicating whether your data are counts.
#' @param family the distribution function to be used in the glm model.
#'   It must be specified as a function: \code{gaussian()}, \code{poisson()}, \code{negative.binomial(theta)}...
#'   If NULL, the family will be \code{negative.binomial(theta)} when \code{counts = TRUE} or \code{gaussian()} when \code{counts = FALSE}.
#' @param theta theta parameter for negative.binomial family.
#' @param epsilon argument to pass to \code{glm.control}, convergence tolerance in the iterative process to estimate the glm model.
#' @param item Name of the analyzed item to show on the screen while \code{sc.p.vector} is in process.
#' @param verbose Name of the analyzed item to show on the screen while \code{T.fit} is in process.
#' @param offset Whether ro use offset for normalization
#'
#' @details \code{rownames(design)} and \code{colnames(data)} must be identical vectors
#'   and indicate array naming. \code{rownames(data)} should contain unique gene IDs.
#'   \code{colnames(design)} are the given names for the variables in the regression model.
#'
#' @return A list containing:
#' \item{SELEC}{matrix containing the expression values for significant genes}
#' \item{sc.p.vector}{vector containing the computed p-values}
#' \item{G}{total number of input genes}
#' \item{g}{number of genes taken in the regression fit}
#' \item{FDR}{p-value at FDR \code{Q} control when Benjamini & Hochberg (BH) correction is used}
#' \item{i}{number of significant genes}
#' \item{dis}{design matrix used in the regression fit}
#' \item{dat}{matrix of expression value data used in the regression fit}
#' \item{...}{additional values from input parameters}
#'
#' @references Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
#' Bioinformatics 22, 1096-1102
#'
#' @author Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}
#'
#' @seealso \code{\link{T.fit}}, \code{\link{lm}}
#'
#' @examples
#' #### GENERATE TIME COURSE DATA
#' ## generates n random gene expression profiles of a data set with
#' ## one control plus 3 treatments, 3 time points and r replicates per time point.
#'
#' # ... (add example usage of the function here)
#'
#' @keywords regression
#'
#' @export
#'
sc.p.vector <- function(scmpObj, Q = 0.05, MT.adjust = "BH", min.obs = 6,
                        counts = FALSE, family = NULL, theta = 10, epsilon = 0.00001,
                        item = "gene", verbose = TRUE, offset = T) {
  # Check the type of the 'design' parameter and set the corresponding variables
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigProClass'"
  )

  # Extract from s4
  dis <- as.data.frame(scmpObj@edesign@dis)
  groups.vector <- scmpObj@edesign@groups.vector
  edesign <- scmpObj@edesign@edesign

  # If 'family' is NULL, determine the distribution function based on 'counts' parameter
  if (is.null(family)) {
    if (!counts) {
      family <- gaussian()
    }
    if (counts) {
      family <- negative.binomial(theta)
    }
  }

  # Convert 'scmpObj' to matrix and select relevant columns based on 'design' rows
  dat <- as.matrix(scmpObj@compress.sce@assays@data@listData$bulk.counts)
  dat <- dat[, as.character(rownames(dis))]
  G <- nrow(dat)

  # Removing rows with many missings:
  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min.obs, ]

  # Removing rows with all zeros:
  sumatot <- apply(dat, 1, sum)
  counts0 <- which(sumatot == 0)
  if (length(counts0) > 0) {
    dat <- dat[-counts0, ]
  }

  # Get dimensions for the input
  g <- dim(dat)[1]
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  sc.p.vector <- vector(mode = "numeric", length = g)

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = g, style = 3)
  }

  # Calculate  offset
  if (offset) {
    offsetData <- log(estimateSizeFactorsForMatrix(dat + 1))
  } else {
    offsetData <- NULL
  }

  # Iterate through each gene and perform the regression fit
  # Store the p-values in 'sc.p.vector'
  for (i in 1:g) {
    y <- as.numeric(dat[i, ])

    # Print progress every 100 genes
    div <- c(1:round(g / 100)) * 100
    if (is.element(i, div)) {
      if (verbose) {
        setTxtProgressBar(pb, i)
      }
      # print(paste(c("fitting ", item, i, "out of", g), collapse = " "))
    }

    model.glm <- glm(y ~ .,
      data = dis, family = family, epsilon = epsilon,
      offset = offsetData
    )
    if (model.glm$null.deviance == 0) {
      sc.p.vector[i] <- 1
    } else {
      model.glm.0 <- glm(y ~ 1,
        family = family, epsilon = epsilon,
        offset = offsetData
      )

      # Perform ANOVA or Chi-square test based on the distribution
      if (family$family == "gaussian") {
        test <- anova(model.glm.0, model.glm, test = "F")
        if (is.na(test[6][2, 1])) {
          sc.p.vector[i] <- 1
        } else {
          sc.p.vector[i] <- test[6][2, 1]
        }
      } else {
        test <- anova(model.glm.0, model.glm, test = "Chisq")
        if (is.na(test[5][2, 1])) {
          sc.p.vector[i] <- 1
        } else {
          sc.p.vector[i] <- test[5][2, 1]
        }
      }
    }
  }

  #----------------------------------------------------------------------
  # Correct p-values using FDR correction and select significant genes
  p.adjusted <- p.adjust(sc.p.vector, method = MT.adjust, n = length(sc.p.vector))
  genes.selected <- rownames(dat)[which(p.adjusted <= Q)]
  FDR <- sort(sc.p.vector)[length(genes.selected)]

  # Subset the expression values of significant genes
  SELEC <- as.matrix(as.data.frame(dat)[genes.selected, ])
  if (nrow(SELEC) == 0) {
    print("no significant genes")
  }

  # Prepare 'sc.p.vector' for output
  sc.p.vector <- as.matrix(sc.p.vector)
  rownames(sc.p.vector) <- rownames(dat)
  colnames(sc.p.vector) <- c("p.value")

  # Add Data to the class
  scPVector.obj <- new("scPVectorClass",
    SELEC = SELEC,
    sc.p.vector = sc.p.vector,
    p.adjusted = p.adjusted,
    G = G,
    g = g,
    FDR = FDR,
    i = nrow(SELEC),
    dis = dis,
    dat = dat,
    min.obs = min.obs,
    Q = Q,
    groups.vector = groups.vector,
    edesign = as.matrix(edesign),
    family = family
  )

  # Update Slot
  scmp.obj@scPVector <- scPVector.obj

  return(scmp.obj)
}
