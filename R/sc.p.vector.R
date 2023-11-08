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
#' @param family the distribution function to be used in the glm model.
#'   It must be specified as a function: \code{gaussian()}, \code{poisson()}, \code{negative.binomial(theta)}...
#'   If NULL, the family will be \code{negative.binomial(theta)} when \code{counts = TRUE} or \code{gaussian()} when \code{counts = FALSE}.
#' @param epsilon argument to pass to \code{glm.control}, convergence tolerance in the iterative process to estimate the glm model.
#' @param verbose Name of the analyzed item to show on the screen while \code{T.fit} is in process.
#' @param offset Whether ro use offset for normalization
#' @param parallel description
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
#' @author Ana Conesa and Maria Jose Nueda, \email{mj.nueda@@ua.es}
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
#' @importFrom stats anova dist glm median na.omit p.adjust
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallelly availableCores
#'
#' @export
#'
sc.p.vector <- function(scmpObj, Q = 0.05, MT.adjust = "BH", min.obs = 6,
                        family = MASS::negative.binomial(theta = 1), epsilon = 0.00001,
                        verbose = TRUE, offset = T, parallel = FALSE) {
  # Check the type of the 'design' parameter and set the corresponding variables
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigProClass'"
  )

  # Extract from s4
  dis <- as.data.frame(scmpObj@edesign@dis)
  groups.vector <- scmpObj@edesign@groups.vector
  edesign <- scmpObj@edesign@edesign

  # Convert 'scmpObj' to matrix and select relevant columns based on 'design' rows
  dat <- as.matrix(scmpObj@compress.sce@assays@data@listData$bulk.counts)
  dat <- dat[, as.character(rownames(dis))]
  G <- nrow(dat)

  # Removing rows with many missings:
  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min.obs, ]

  # Add check
  # assert_that((dat@Dim[1] > 1), msg = paste(min.obs, "for 'min.obs' is too high. Try lowering the threshold."))
  assert_that((nrow(dat) > 1), msg = paste(min.obs, "for 'min.obs' is too high. Try lowering the threshold."))

  # if(verbose){
  #     message(paste("'min.obs' is set at", min.obs))
  #     message(paste("After filtering with 'min.obs'", scmpObj@compress.sce@assays@data@listData$bulk.counts@Dim[1] - dat@Dim[1], "gene are dropped"))
  # }

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

  if (parallel == F) {
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = g, style = 3)
    }
  }

  # Calculate  offset
  if (offset) {
    dat <- dat + 1
    offsetData <- log(scmp_estimateSizeFactorsForMatrix(dat))
    if (verbose) {
      message("Using DESeq2::estimateSizeFactorsForMatrix")
      message("Please cite DESeq2 as 'Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)'")
    }
  } else {
    offsetData <- NULL
  }

  if (parallel) {
    numCores <- availableCores()
    if (verbose) {
      message(paste("Running with", numCores, "cores..."))
    }
  } else {
    numCores <- 1
  }

  p.vector.list <- mclapply(1:g, function(i, g_lapply = g, dat_lapply = dat, dis_lapply = dis, family_lapply = family, epsilon_lapply = epsilon, offsetdata_lapply = offsetData, pb_lapply = pb, verbose_lapply = verbose) {
    y <- as.numeric(dat_lapply[i, ])

    # Print prog_lapplyress every 100 g_lapplyenes
    div <- c(1:round(g_lapply / 100)) * 100

    if (parallel == F) {
      if (is.element(i, div) && verbose_lapply) {
        if (verbose) {
          setTxtProgressBar(pb_lapply, i)
        }
      }
    }
    model.glm <- glm(y ~ .,
      data = dis_lapply, family = family_lapply, epsilon = epsilon_lapply,
      offset = offsetdata_lapply
    )

    sc_p_val <- NA
    if (model.glm$null.deviance == 0) {
      sc_p_val <- 1
    } else {
      model.glm.0 <- glm(y ~ 1,
        family = family_lapply, epsilon = epsilon_lapply,
        offset = offsetdata_lapply
      )

      # Perform ANOVA or Chi-square test based on the dis_lapplytribution
      if (family_lapply$family == "gaussian") {
        test <- anova(model.glm.0, model.glm, test = "F")
        sc_p_val <- ifelse(is.na(test[6][2, 1]), 1, test[6][2, 1])
      } else {
        test <- anova(model.glm.0, model.glm, test = "Chisq")
        sc_p_val <- ifelse(is.na(test[5][2, 1]), 1, test[5][2, 1])
      }
    }

    return(sc_p_val)
  }, mc.cores = numCores)

  names(p.vector.list) <- rownames(dat)
  sc.p.vector <- unlist(p.vector.list, recursive = T, use.names = T)
  #----------------------------------------------------------------------
  # Correct p-values using FDR correction and select significant genes
  p.adjusted <- unlist(p.adjust(sc.p.vector, method = MT.adjust, n = length(sc.p.vector)),
    recursive = T, use.names = T
  )
  names(p.adjusted) <- names(sc.p.vector)
  genes.selected <- rownames(dat)[which(p.adjusted <= Q)]
  FDR <- sort(sc.p.vector)[length(genes.selected)]

  # Subset the expression values of significant genes
  SELEC <- dat[rownames(dat) %in% genes.selected, ]

  if (nrow(SELEC) == 0) {
    message("No significant genes detected. Try changing parameters.")
    return(scmpObj)
  } else {
    # Prepare 'sc.p.vector' for output
    names(sc.p.vector) <- rownames(dat)

    # Add Data to the class
    scPVector.obj <- new("scPVectorClass",
      SELEC = as(SELEC, "dgCMatrix"),
      p.vector = sc.p.vector,
      p.adjusted = p.adjusted,
      FDR = FDR,
      dis = dis,
      groups.vector = groups.vector,
      family = family
    )

    # Update Slot
    scmpObj@scPVector <- scPVector.obj

    # Update Parameter Slot
    scmpObj@addParams@Q <- Q
    scmpObj@addParams@min.obs <- min.obs
    scmpObj@addParams@g <- g
    scmpObj@addParams@MT.adjust <- MT.adjust
    scmpObj@addParams@epsilon <- epsilon
    scmpObj@distribution <- family

    return(scmpObj)
  }
}
