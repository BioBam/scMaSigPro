#' @title Perform fitting with full model.
#'
#' @description
#' Performs a regression fit for each gene taking all variables present in the
#' model.
#'
#' @importFrom stats anova dist glm median na.omit p.adjust glm.control
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallelly availableCores
#' @importFrom MASS negative.binomial glm.nb
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param p_value Significance level used for variable selection in the stepwise
#' regression.
#' @param mt_correction A character string specifying the p-value correction
#' method.
#' @param min_na Minimum values needed per gene across cells to estimate the
#' model.
#' @param family  Distribution of the error term.
#' @param link Type of link function to use in the model. Default is "log".
#' @param epsilon Model convergence tolerance.
#' @param offset logical value specifying whether to use offset during fitting.
#' @param log_offset A logical value specifying whether to take the logarithm of
#' the offsets.
#' @param max_it Maximum number of iterations to fit the model.
#' @param Type of link function to use in the model. Default is "log".
#' @param parallel Use forking process to run parallelly. (Default is FALSE)
#' (Currently, Windows is not supported)
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with updated `Profile`
#' slot.
#'
#' @seealso \code{\link{VariableProfiles}}
#'
#' @references Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles
#' in Time-Course Microarray Experiments. Bioinformatics 22, 1096-1102
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#'
#' @keywords regression models
#' @export
sc.p.vector <- function(scmpObj, p_value = 0.05, mt_correction = "BH",
                        min_na = 6,
                        family = negative.binomial(theta = 10),
                        epsilon = 1e-8,
                        verbose = TRUE,
                        offset = TRUE,
                        parallel = FALSE,
                        log_offset = FALSE,
                        max_it = 100,
                        link = "log") {
  # Check the type of the 'design' parameter and set the corresponding variables
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'ScMaSigPro'"
  )

  # Extract from s4
  dis <- as.data.frame(scmpObj@Design@predictor_matrix)
  groups.vector <- scmpObj@Design@groups.vector
  alloc <- scmpObj@Design@assignment_matrix

  # Convert 'scmpObj' to matrix and select relevant columns based on 'design' rows
  dat <- as.matrix(scmpObj@Dense@assays@data@listData$bulk.counts)
  dat <- dat[, as.character(rownames(dis))]
  G <- nrow(dat)

  # Check for the log function
  assert_that(link %in% c("log", "identity"),
    msg = "link function should be either 'log' or 'identity'"
  )

  # Update the family
  family[["link"]] <- link

  # Add check
  # assert_that((dat@Dim[1] > 1), msg = paste(min_na, "for 'min_na' is too high. Try lowering the threshold."))
  assert_that(min_na <= ncol(dat),
    msg = paste(
      min_na,
      "for 'min_na' is too high. Try lowering the threshold."
    )
  )

  # Removing rows with many missings:
  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min_na, ]
  # if(verbose){
  #     message(paste("'min_na' is set at", min_na))
  #     message(paste("After filtering with 'min_na'", scmpObj@Dense@assays@data@listData$bulk.counts@Dim[1] - dat@Dim[1], "gene are dropped"))
  # }

  # Removing rows with all zeros:
  # dat[is.na(dat)] <- 0
  # dat <- dat[rowSums(dat) != 0, , drop = FALSE]
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

  if (parallel == FALSE) {
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = g, style = 3)
    }
  }

  # Calculate  offset
  if (offset) {
    dat <- dat + 1
    offsetData <- scmp_estimateSizeFactorsForMatrix(dat)
    if (log_offset) {
      offsetData <- log(offsetData)
    }
    scmpObj@Design@offset <- offsetData
  } else {
    offsetData <- NULL # scmpObj@Design@offset
  }

  if (parallel) {
    os_name <- get_os()
    if (os_name == "windows") {
      numCores <- 1
      warning("Currently, we only support sequential processing on windows based systems...")
    } else {
      numCores <- availableCores() - 1
    }
    if (verbose) {
      message(paste("Running with", numCores, "cores..."))
    }
  } else {
    numCores <- 1
  }

  # Check for weight usage
  useWeights <- FALSE
  logWeights <- FALSE
  useInverseWeights <- FALSE
  if (useWeights) {
    # Get the pathframe
    compressed.data <- as.data.frame(scmpObj@Dense@colData)

    # Get bin_name and bin size
    weight_df <- compressed.data[, c(scmpObj@Parameters@bin_size_colname),
      drop = TRUE
    ]

    # Set names
    names(weight_df) <- rownames(compressed.data)

    # Log of weights
    if (logWeights) {
      weight_df <- log(weight_df)
    }

    # Take inverse weights
    if (useInverseWeights) {
      weight_df <- 1 / weight_df
    }
  } else {
    weight_df <- NULL
  }
  p.vector.list <- mclapply(1:g, function(i, g_lapply = g, dat_lapply = dat,
                                          dis_lapply = dis,
                                          family_lapply = family,
                                          epsilon_lapply = epsilon,
                                          offsetdata_lapply = offsetData,
                                          pb_lapply = pb,
                                          weights_lapply = weight_df,
                                          verbose_lapply = verbose,
                                          max_it_lapply = max_it) {
    y <- as.numeric(dat_lapply[i, ])

    # Print prog_lapplyress every 100 g_lapplyenes
    div <- c(1:round(g_lapply / 100)) * 100

    if (parallel == FALSE) {
      if (is.element(i, div) && verbose_lapply) {
        if (verbose) {
          setTxtProgressBar(pb_lapply, i)
        }
      }
    }

    # Set full Model
    model.glm <- glm(y ~ .,
      data = dis_lapply, family = family_lapply, epsilon = epsilon_lapply,
      offset = offsetdata_lapply, weights = weights_lapply,
      maxit = max_it_lapply
    )
    sc_p_val <- NA
    if (model.glm$null.deviance == 0) {
      sc_p_val <- 1
    } else {
      model.glm.0 <- glm(y ~ 1,
        family = family_lapply, epsilon = epsilon_lapply,
        offset = offsetdata_lapply, weights = weights_lapply,
        maxit = max_it_lapply
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
  # Stop the cluster

  names(p.vector.list) <- rownames(dat)
  sc.p.vector <- unlist(p.vector.list, recursive = T, use.names = T)
  #----------------------------------------------------------------------
  # Correct p-values using FDR correction and select significant genes
  p.adjusted <- unlist(
    p.adjust(sc.p.vector,
      method = mt_correction,
      n = length(sc.p.vector)
    ),
    recursive = T, use.names = T
  )
  names(p.adjusted) <- names(sc.p.vector)
  genes.selected <- rownames(dat)[which(p.adjusted <= p_value)]
  FDR <- sort(sc.p.vector)[length(genes.selected)]

  # Subset the expression values of significant genes
  SELEC <- dat[rownames(dat) %in% genes.selected, , drop = FALSE]

  if (nrow(SELEC) == 0) {
    message("No significant genes detected. Try changing parameters.")
    return(scmpObj)
  } else {
    # Prepare 'sc.p.vector' for output
    names(sc.p.vector) <- rownames(dat)

    # Add Data to the class
    profile.obj <- new("VariableProfiles",
      non_flat = rownames(SELEC),
      p_values = sc.p.vector,
      adj_p_values = p.adjusted,
      fdr = FDR
    )

    # Update Slot
    scmpObj@Profile <- profile.obj

    # Update Parameter Slot useInverseWeights
    # scmpObj@Parameters@useWeights <- useWeights
    scmpObj@Parameters@log_offset <- log_offset
    # scmpObj@Parameters@logWeights <- logWeights
    scmpObj@Parameters@max_it <- as.integer(max_it)
    # scmpObj@Parameters@useInverseWeights <- useInverseWeights
    scmpObj@Parameters@offset <- offset
    scmpObj@Parameters@p_value <- p_value
    scmpObj@Parameters@min_na <- min_na
    scmpObj@Parameters@g <- g
    scmpObj@Parameters@mt_correction <- mt_correction
    scmpObj@Parameters@epsilon <- epsilon
    scmpObj@Parameters@distribution <- family
    scmpObj@Parameters@link <- link

    return(scmpObj)
  }
}
