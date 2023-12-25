#' Makes a stepwise regression fit for time series gene expression experiments
#'
#' \code{s.t.fit} selects the best regression model for each gene using stepwise regression.
#'
#' @param scmpObj Can either be a \code{\link{p.vector}} object or a matrix containing expression scmpObj with the same requirements as for
#' the \code{\link{p.vector}} function.
#' @param selection_method Argument to be passed to the step function. Can be either \code{"backward"}, \code{"forward"}, \code{"two.ways.backward"}, or \code{"two.ways.forward"}.
#' @param p_value Significance level used for variable selection in the stepwise regression.
#' @param nvar_correction Argument for correcting T.fit significance level. See details.
#' @param family The distribution function to be used in the glm model. It must be the same used in \code{p.vector}.
#' @param epsilon Argument to pass to \code{glm.control}, convergence tolerance in the iterative process to estimate the glm model.
#' @param verbose Name of the analyzed item to show on the screen while \code{T.fit} is in process.
#' @param offset Whether ro use offset for normalization
#' @param parallel description
#' @param log_offset description
#' @param max_it description
#'
#' @details
#' In the maSigPro approach, \code{\link{p.vector}} and \code{\link{T.fit}} are subsequent steps, meaning that significant genes are
#' first selected based on a general model, and then the significant variables for each gene are found by step-wise regression.
#'
#' The step regression can be \code{"backward"} or \code{"forward"}, indicating whether the step procedure starts from the
#' model with all or none variables. With the \code{"two.ways.backward"} or \code{"two.ways.forward"} options, the variables are both allowed to get in and out.
#' At each step, the p-value of each variable is computed, and variables get in/out of the model when this p-value is
#' lower or higher than the given threshold \code{p_value}. When \code{nvar_correction} is TRUE, the given significance level is corrected by the number of variables in the model.
#'
#' @return
#' A list containing the following elements:
#' \item{sol}{Matrix for summary results of the stepwise regression. For each selected gene, the following values are given:
#' \itemize{
#'   \item p-value of the regression ANOVA
#'   \item R-squared of the model
#'   \item p-value of the regression coefficients of the selected variables
#' }}
#' \item{coefficients}{Matrix containing regression coefficients for the adjusted models.}
#' \item{group.coeffs}{Matrix containing the coefficients of the implicit models of each experimental group.}
#' \item{variables}{Variables in the complete regression model.}
#' \item{G}{Total number of input genes.}
#' \item{g}{Number of genes taken in the regression fit.}
#' \item{dat}{Input analysis scmpObj matrix.}
#' \item{dis}{Regression design matrix.}
#' \item{selection_method}{Imputed step method for stepwise regression.}
#' \item{alloc}{Matrix of experimental design.}
#' \item{influ.info}{scmpObj frame of genes containing influential scmpObj.}
#'
#' @references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
#' Bioinformatics 22, 1096-1102}
#'
#' @author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@@ua.es}}
#'
#' @seealso{\code{\link{p.vector}}, \code{\link{step}}}
#'
#' @importFrom maSigPro position reg.coeffs
#' @importFrom stats influence.measures
#' @keywords regression
#' @keywords models
#' @export
sc.t.fit <- function(scmpObj,
                     selection_method = "backward",
                     p_value = scmpObj@Parameters@p_value,
                     nvar_correction = FALSE,
                     family = scmpObj@Parameters@distribution,
                     epsilon = scmpObj@Parameters@epsilon,
                     offset = scmpObj@Parameters@offset,
                     verbose = TRUE,
                     parallel = FALSE,
                     log_offset = scmpObj@Parameters@log_offset,
                     max_it = scmpObj@Parameters@max_it) {
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'ScMaSigPro'"
  )

  # Transfer Data
  dis <- scmpObj@Design@predictor_matrix
  p_value <- scmpObj@Parameters@p_value
  groups.vector <- scmpObj@Design@groups.vector
  groups.vector <- c(groups.vector[nchar(groups.vector) == min(nchar(groups.vector))][1], groups.vector)
  alloc <- scmpObj@Design@assignment_matrix
  G <- scmpObj@Parameters@g

  dat.all <- scmpObj@Dense@assays@data@listData$bulk.counts
  dat <- dat.all[scmpObj@Profile@non_flat, , drop = F]
  dat <- rbind(c(rep(1, ncol(dat))), dat)
  dat <- dat[, as.character(rownames(dis))]

  g <- (dim(dat)[1] - 1)
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  vars.in <- colnames(dis)
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
  rownames(influ.info) <- rownames(dis)

  if (nvar_correction) {
    p_value <- p_value / ncol(dis)
  }

  # Check for weight usage
  useWeights <- FALSE
  useInverseWeights <- FALSE
  logWeights <- FALSE
  if (useWeights) {
    # Get the pathframe
    compressed.data <- as.data.frame(scmpObj@Dense@colData)

    # Get bin_name and bin size
    weight_df <- compressed.data[, c(scmpObj@Parameters@bin_size_colname), drop = TRUE]

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

  # Calculate  offset
  if (offset) {
    dat <- dat + 1
    offsetData <- scmp_estimateSizeFactorsForMatrix(dat)
    if (log_offset) {
      offsetData <- log(offsetData)
    }
  } else {
    offsetData <- NULL
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
    if (verbose) {
      i <- 1
      div <- c(1:round(g / 100)) * 100
      pb <- txtProgressBar(min = 0, max = g, style = 3)
    }
  }

  # Collect all the y
  y_input <- parallel::mclapply(2:(g + 1), function(i, dat_lapply = dat) {
    return(as.numeric(dat_lapply[i, ]))
  }, mc.cores = numCores)
  names(y_input) <- rownames(dat)[-1]

  # Select the covariates
  if (selection_method == "backward") {
    result_list <- parallel::mclapply(names(y_input), function(gene_name, dat_lapply = dat, dis_lapply = dis, family_lapply = family, epsilon_lapply = epsilon, offsetData_lapply = offsetData, pb_lapply = pb, verbose_lapply = verbose, vars_in_lapply = vars.in, Q_lapply = p_value, influ.info_lapply = influ.info, weights_lapply = weight_df, max_it_lapply = max_it) {
      # result_list <- lapply(names(y_input), function(gene_name, g_lapply = g, dat_lapply = dat, dis_lapply = dis, family_lapply = family, epsilon_lapply = epsilon, offsetData_lapply = offsetData, pb_lapply = pb, verbose_lapply = verbose, vars_in_lapply = vars.in, Q_lapply = Q, influ.info_lapply = influ.info) {
      y <- y_input[[gene_name]]


      reg_scmpObj <- sc.stepback(y = y, d = as.data.frame(dis_lapply), alfa = Q_lapply, family = family_lapply, epsilon = epsilon_lapply, useOffset = offsetData_lapply, useWeight = weights_lapply, max_it = max_it_lapply)
      lmf_scmpObj <- glm(y ~ ., data = as.data.frame(dis_lapply), family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      model.glm.0_scmpObj <- glm(y ~ 1, family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      if (parallel == FALSE) {
        if (verbose_lapply) {
          if (verbose_lapply) {
            i <- i + 1
            setTxtProgressBar(pb_lapply, i)
          }
        }
      }
      return(extract_fitting(reg = reg_scmpObj, lmf = lmf_scmpObj, model.glm.0 = model.glm.0_scmpObj, dis = dis_lapply, family = family_lapply, name = gene_name, vars.in = vars_in_lapply, alfa = Q_lapply, influ.info = influ.info_lapply))
      # return(list(
      #     reg = reg,
      #     lmf = lmf,
      #     model.glm.0 = model.glm.0
      # ))
    }, mc.cores = numCores, mc.set.seed = 2023)
    # })
  } else if (selection_method == "forward") {
    result_list <- parallel::mclapply(names(y_input), function(gene_name, g_lapply = g, dat_lapply = dat, dis_lapply = dis, family_lapply = family, epsilon_lapply = epsilon, offsetData_lapply = offsetData, pb_lapply = pb, verbose_lapply = verbose, vars_in_lapply = vars.in, Q_lapply = p_value, influ.info_lapply = influ.info, weights_lapply = weight_df, max_it_lapply = max_it) {
      y <- y_input[[gene_name]]

      reg_scmpObj <- sc.stepfor(y = y, d = as.data.frame(dis_lapply), alfa = Q_lapply, family = family_lapply, epsilon = epsilon_lapply, useOffset = offsetData_lapply, useWeight = weights_lapply, max_it = max_it_lapply)
      lmf_scmpObj <- glm(y ~ ., data = as.data.frame(dis_lapply), family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      model.glm.0_scmpObj <- glm(y ~ 1, family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      div <- c(1:round(g / 100)) * 100
      if (parallel == FALSE) {
        if (is.element(y, div) && verbose_lapply) {
          if (verbose) {
            setTxtProgressBar(pb_lapply, y)
          }
        }
      }
      return(extract_fitting(reg = reg_scmpObj, lmf = lmf_scmpObj, model.glm.0 = model.glm.0_scmpObj, dis = dis_lapply, family = family_lapply, name = gene_name, vars.in = vars_in_lapply, alfa = Q_lapply, influ.info = influ.info_lapply))
    }, mc.cores = numCores, mc.set.seed = 2023)
  } else if (selection_method == "two.ways.backward") {
    result_list <- parallel::mclapply(names(y_input), function(gene_name, g_lapply = g, dat_lapply = dat, dis_lapply = dis, family_lapply = family, epsilon_lapply = epsilon, offsetData_lapply = offsetData, pb_lapply = pb, verbose_lapply = verbose, vars_in_lapply = vars.in, Q_lapply = p_value, influ.info_lapply = influ.info, weights_lapply = weight_df, max_it_lapply = max_it) {
      y <- y_input[[gene_name]]

      reg_scmpObj <- sc.two.ways.stepback(y = y, d = as.data.frame(dis_lapply), alfa = Q_lapply, family = family_lapply, epsilon = epsilon_lapply, useOffset = offsetData_lapply, useWeight = weights_lapply, max_it = max_it_lapply)
      lmf_scmpObj <- glm(y ~ ., data = as.data.frame(dis_lapply), family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      model.glm.0_scmpObj <- glm(y ~ 1, family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      div <- c(1:round(g / 100)) * 100
      if (parallel == FALSE) {
        if (is.element(y, div) && verbose_lapply) {
          if (verbose) {
            setTxtProgressBar(pb_lapply, y)
          }
        }
      }
      return(extract_fitting(reg = reg_scmpObj, lmf = lmf_scmpObj, model.glm.0 = model.glm.0_scmpObj, dis = dis_lapply, family = family_lapply, name = gene_name, vars.in = vars_in_lapply, alfa = Q_lapply, influ.info = influ.info_lapply))
    }, mc.cores = numCores, mc.set.seed = 2023)
    # })
  } else if (selection_method == "two.ways.forward") {
    result_list <- parallel::mclapply(names(y_input), function(gene_name, g_lapply = g, dat_lapply = dat, dis_lapply = dis, family_lapply = family, epsilon_lapply = epsilon, offsetData_lapply = offsetData, pb_lapply = pb, verbose_lapply = verbose, vars_in_lapply = vars.in, Q_lapply = p_value, influ.info_lapply = influ.info, weights_lapply = weight_df, max_it_lapply = max_it) {
      y <- y_input[[gene_name]]

      reg_scmpObj <- sc.two.ways.stepfor(y = y, d = as.data.frame(dis_lapply), alfa = Q_lapply, family = family_lapply, epsilon = epsilon_lapply, useOffset = offsetData_lapply, useWeight = weights_lapply, max_it = max_it_lapply)
      lmf_scmpObj <- glm(y ~ ., data = as.data.frame(dis_lapply), family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      model.glm.0_scmpObj <- glm(y ~ 1, family = family_lapply, epsilon = epsilon_lapply, offset = offsetData_lapply, weights = weights_lapply, maxit = max_it_lapply)
      div <- c(1:round(g / 100)) * 100
      if (parallel == FALSE) {
        if (is.element(y, div) && verbose_lapply) {
          if (verbose) {
            setTxtProgressBar(pb_lapply, y)
          }
        }
      }
      return(extract_fitting(reg = reg_scmpObj, lmf = lmf_scmpObj, model.glm.0 = model.glm.0_scmpObj, dis = dis_lapply, family = family_lapply, name = gene_name, vars.in = vars_in_lapply, alfa = Q_lapply, influ.info = influ.info_lapply))
    }, mc.cores = numCores, mc.set.seed = 2023)
  } else {
    stop("stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward")
  }

  #--------------------------------------------
  feature_names <- unlist(lapply(result_list, function(element) {
    return(element[["feature_name"]])
  }))
  # Get the soluction frame
  sol.list <- lapply(result_list, function(element) {
    return(element[["sol"]])
  })
  # Get Coeffcient
  coeff.list <- lapply(result_list, function(element) {
    return(element[["coeff"]])
  })
  # Get t scores
  t.list <- lapply(result_list, function(element) {
    return(element[["t"]])
  })
  # Get influ.info
  influ.info.list <- lapply(result_list, function(element) {
    return(element[["influ.info"]])
  })

  # Assuming 'parallel' is your list
  # influ.info.list <- influ.info.list[!sapply(influ.info.list, function(x) is.logical(x))]
  influ.info.list <- influ.info.list[!vapply(influ.info.list, is.logical, logical(1))]

  # Lapply to remove column 1
  influ.info.list <- lapply(influ.info.list, function(element) {
    return(element[, -1, drop = FALSE])
  })

  # Create scmpObjframe
  sol <- do.call("rbind", sol.list)
  coefficients <- do.call("rbind", coeff.list)
  t.score <- do.call("rbind", t.list)
  influ.info <- do.call("cbind", influ.info.list)

  # Add rownames
  rownames(coefficients) <- feature_names
  rownames(t.score) <- feature_names

  #-------------------------
  # Ends here

  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    coeffic <- coefficients
    t.score <- as.data.frame(t.score)
    colnames(sol) <- c(
      "p-value", "R-squared", "p.valor_beta0",
      paste("p.valor_", vars.in, sep = "")
    )
    colnames(coefficients) <- c("beta0", paste("beta", vars.in,
      sep = ""
    ))
    colnames(t.score) <- c("t.score_beta0", paste("t.score_",
      vars.in,
      sep = ""
    ))
    if (!is.null(groups.vector) & !is.null(alloc)) {
      groups <- colnames(alloc)[3:ncol(alloc)]
      degree <- (length(groups.vector) / length(groups)) -
        1
      for (w in 1:nrow(coefficients)) {
        A <- NULL
        col.names <- NULL
        for (l in 1:length(groups)) {
          B <- reg.coeffs(coefficients = coefficients[w, ], groups.vector = groups.vector, group = groups[l])
          cols <- paste(rep(groups[l], each = length(B)),
            paste("beta", c(0:(length(B) - 1)), sep = ""),
            sep = "_"
          )
          A <- c(A, B)
          col.names <- c(col.names, cols)
        }
        group.coeffs <- (rbind(group.coeffs, A))
      }
      colnames(group.coeffs) <- col.names
      rownames(group.coeffs) <- rownames(coefficients)
    }
  }

  if (!is.null(influ.info)) {
    if (verbose) {
      message(paste("\nInfluence:", ncol(influ.info), "genes with influential onservations detected. Model validation for these genes is recommended"))
    }
  } else {
    if (verbose) {
      message(paste("\nNo genes with influential observations detected"))
    }
    influ.info <- matrix(data = NA, nrow = 0, ncol = 0)
  }
  # influ.info <- influ.info[, -1]

  # Create a constructor for the class
  t.fit.object <- new("Estimates",
    significance_matrix = as.matrix(sol),
    coefficient_matrix = as.matrix(coefficients),
    path_coefficient_matrix = group.coeffs,
    t_score_matrix = as.matrix(t.score),
    path = groups.vector,
    influential = influ.info
  )

  # Added Tfit
  scmpObj@Estimate <- t.fit.object

  # Update Parameter Slot
  scmpObj@Parameters@p_value <- p_value
  scmpObj@Parameters@log_offset <- log_offset
  scmpObj@Parameters@max_it <- as.integer(max_it)
  scmpObj@Parameters@offset <- offset
  scmpObj@Parameters@epsilon <- epsilon
  scmpObj@Parameters@selection_method <- selection_method
  scmpObj@Parameters@distribution <- family

  return(scmpObj)
}
