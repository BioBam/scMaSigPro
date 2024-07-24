## Show functions for scMaSigPro
## Description: Following functions are used to extract information from the ScMaSigProClass
## Object
# 1. showCoeff():
# 2. showInflu():
# 3. showTS():
# 4. showSol():
# 5. showSigProf():
# 6. .ScMaSigPro_show(): For object cat
# 7. extract_info():
# 8. showParams():
# 9. showPoly():
# 10. showGroupCoeff():

###############################################################################
#' @title Show or Return the Coefficient matrix
#'
#' @description
#' This function is used to view or return the Coefficients from the provided
#' scMaSigPro object.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param includeInflu Whether to include genes with influential observations.
#' @param return Whether to return the data. (Default: TRUE)
#'
#' @return The computed Coefficient matrix as a dataframe.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showCoeff <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@Estimate@coefficient_matrix) == c(0, 0)),
    msg = "'coefficient_matrix' is not computed yet"
  )

  # Extract
  coefficients <- scmpObj@Estimate@coefficient_matrix %>% as.data.frame()

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    coefficients <- coefficients[!(rownames(coefficients) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(coefficients)
  }

  # If requested
  if (return) {
    return(coefficients)
  }
}

###############################################################################
#' @title Return the matrix of genes with influential observation
#'
#' @description
#' This function is used to view or return the matrix of genes with influential
#' observation from the provided scMaSigPro object.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param return Whether to return the data. (Default: TRUE)
#'
#' @return Matrix of genes with influential observation.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showInflu <- function(scmpObj, view = FALSE, return = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@Estimate@influential) == c(0, 0)),
    msg = "tscore is not computed yet"
  )

  # Extract
  influ <- scmpObj@Estimate@influential

  # If viewing is requested
  if (view) {
    View(influ)
  }

  # If requested
  if (return) {
    return(influ)
  }
}

###############################################################################
#' @title Show or Return the t-score matrix
#'
#' @description
#' This function is used to view or return the t-scores from the provided
#' scMaSigPro object.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param includeInflu Whether to include genes with influential observations.
#' @param return Whether to return the data. (Default: TRUE)
#'
#' @return The computed t-score matrix as a dataframe.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showTS <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@Estimate@t_score_matrix) == c(0, 0)),
    msg = "tscore is not computed yet"
  )

  # Extract
  tscore <- scmpObj@Estimate@t_score_matrix

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    tscore <- tscore[!(rownames(tscore) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(tscore)
  }

  # If requested
  if (return) {
    return(tscore)
  }
}

###############################################################################
#' @title Show or Return the P-values after model fitting.
#'
#' @description
#' This function is used to view or return the matrix of p-values for each term
#' and the full model from the provided scMaSigPro object.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param includeInflu Whether to include genes with influential observations.
#' @param return Whether to return the data. (Default: TRUE)
#'
#' @return The computed p-values for each term and full model as a dataframe.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showSol <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@Estimate@significance_matrix) == c(0, 0)),
    msg = "'significance_matrix' is not computed yet"
  )

  # Extract
  sol <- scmpObj@Estimate@significance_matrix %>% as.data.frame()

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    sol <- sol[!(rownames(sol) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(sol)
  }

  # If requested
  if (return) {
    return(sol)
  }
}

###############################################################################
#' @title Show or Return the counts for non-flat profile.
#'
#' @description
#' This function is used to view or return the pseudo-bulk counts of the genes
#' with non-flat profiles. from the provided scMaSigPro object.
#'
#' @importFrom utils View
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param includeInflu Whether to include genes with influential observations.
#' @param return Whether to return the data. (Default: TRUE)
#'
#' @return Pseudo-bulk counts as matrix for genes with non-flat profiles.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showSigProf <- function(scmpObj, view = FALSE, return = TRUE,
                        includeInflu = FALSE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@Estimate@significance_matrix) == c(0, 0)),
    msg = "Sol is not computed yet"
  )

  # Extract
  sol <- showSol(scmpObj,
    view = FALSE, return = TRUE,
    includeInflu = includeInflu
  ) %>% as.data.frame()
  # Extract rownames
  bulk.counts <- scmpObj@Dense@assays@data@listData$bulk.counts
  bulk.counts <- bulk.counts[rownames(bulk.counts) %in% rownames(sol), , drop = FALSE]

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    bulk.counts <- bulk.counts[!(rownames(bulk.counts) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(as.matrix(bulk.counts))
  }

  # If requested
  if (return) {
    return(bulk.counts)
  }
}

###############################################################################
#' @title Print the full model formula.
#'
#' @description
#' Print the full model formula in console as a string.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#'
#' @return Character string of the formula for the full model.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showPoly <- function(scmpObj) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(
    all(!is.na(colnames(scmpObj@Design@predictor_matrix)) | length(
      colnames(scmpObj@Design@predictor_matrix) > 1
    )),
    msg = "Please setup the model first, using 'sc.make.design.matrix()'"
  )

  # Extract columns
  df.col <- colnames(scmpObj@Design@predictor_matrix)

  # Extract betas
  beta_names <- paste0("beta", seq(1:length(df.col)))

  # Generate formula string
  formula_parts <- vapply(seq_along(df.col),
    function(i) paste(beta_names[i], "*", df.col[i], sep = ""),
    FUN.VALUE = character(1)
  )

  # Make formula
  formula_string <- paste("beta0", paste(formula_parts, collapse = " + "),
    sep = " + "
  )

  return(formula_string)
}

###############################################################################
#' @title Show the parameters used during the workflow.
#'
#' @description
#' Get or View all the parameters used during the workflow.
#'
#' @importFrom methods slot slotNames
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param return Whether to return the data. (Default: TRUE)
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @return Dataframe of the parameters used in the analysis.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
showParams <- function(scmpObj, view = FALSE, return = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Get slot names, assuming 'param' is a slot within 'scmpObj'
  all_slots <- slotNames(scmpObj)

  # Get 'param' slot data using the correct S4 accessor method
  paramData <- slot(scmpObj, "Parameters")

  # Get all slots of 'param', assuming 'param' itself is an S4 object with slots
  params.frame.list <- lapply(slotNames(paramData), function(parameter) {
    slot(paramData, parameter)
  })
  names(params.frame.list) <- slotNames(paramData)
  params.frame.list[["distribution"]] <- NULL

  # Get the data
  params.frame <- data.frame(
    parameters = names(params.frame.list),
    value = unlist(params.frame.list)
  )

  # Add family
  distribution.f <- data.frame(
    parameters = "distribution",
    value = paste0(scmpObj@Parameters@distribution[["family"]])
  )


  # Add
  params <- rbind(
    params.frame,
    distribution.f
  )
  rownames(params) <- NULL

  # If viewing is requested and the 'View' function is available
  if (view && exists("View")) {
    View(params)
  }

  # If requested, return the parameters
  if (return) {
    return(params)
  }
}

###############################################################################
#' @title Show or Return the Branching Path Coefficient matrix
#'
#' @description
#' This function is used to view or return the branching paths Coefficients from
#' the provided scMaSigPro object.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param view Whether to view the data in the explorer. (Default: FALSE)
#' @param return Whether to return the data. (Default: TRUE)
#' @param includeInflu Whether to include genes with influential observations.
#'
#' @return The computed branching path Coefficient matrix as a dataframe.
#'
#' @export
showGroupCoeff <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@Estimate@path_coefficient_matrix) == c(0, 0)),
    msg = "path_coefficient_matrix is not computed yet"
  )

  # Extract
  grpCoeff <- scmpObj@Estimate@path_coefficient_matrix %>% as.data.frame()

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    grpCoeff <- grpCoeff[!(rownames(grpCoeff) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(grpCoeff)
  }

  # If requested
  if (return) {
    return(grpCoeff)
  }
}

###############################################################################
#' @title Show ScMaSigPro Object Information
#'
#' @description
#' This method displays basic information about the ScMaSigPro object when the
#' object is printed in the console. The method is automatically called when the
#' user writes the name of the object in the console.
#'
#' @param object An object of class \code{ScMaSigPro}.
#'
#' @importFrom S4Vectors coolcat
#'
#' @keywords internal
#' @export
.ScMaSigPro_show <- function(object) {
  # Show Basic information
  cat("Class: ScMaSigProClass\n")
  cat(paste0("nCells: ", ncol(object@Sparse), "\n"))
  cat(paste0("nFeatures: ", nrow(object@Sparse), "\n"))
  cat("Pseudotime Range:", paste(round(
    range(SingleCellExperiment::colData(object@Sparse)[[object@Parameters@ptime_col]]), 3
  )))
  cat(paste("\nBranching Paths:", paste(
    unique(SingleCellExperiment::colData(object@Sparse)[[object@Parameters@path_col]]),
    collapse = ", "
  )))

  # Calculate the Compression
  compressed.cell.metadata <- cDense(object)
  if (length(compressed.cell.metadata) > 0) {
    cat(paste0(
      "\nBinned Pseudotime: ", paste(
        range(
          compressed.cell.metadata[[object@Parameters@bin_ptime_col]]
        ),
        collapse = "-"
      ), "(Range), ",
      round(
        mean(
          compressed.cell.metadata[[object@Parameters@bin_ptime_col]]
        ),
        2
      ), "(Mean), "
    ))

    # Extract info
    per_path_num_bin <- extract_info(compressed.cell.metadata,
      return_type = "num_bins",
      bin_size_col = object@Parameters@bin_size_col,
      object@Parameters@path_col
    )
    per_path_bin_size <- round(extract_info(
      compressed.cell.metadata,
      return_type = "avg_bin_size",
      bin_size_col = object@Parameters@bin_size_col,
      object@Parameters@path_col
    ))

    # Paste
    cat("\nNumber of bins->", paste(names(per_path_num_bin),
      per_path_num_bin,
      sep = ": "
    ))
    cat("\nAverage bin Size->", paste(names(per_path_bin_size),
      per_path_bin_size,
      sep = ": "
    ))
  }

  # Influential Genes if any
  if (length(colnames(object@Design@predictor_matrix)) > 0) {
    cat(paste("\nPolynomial Order:", object@Parameters@poly_degree))
  }

  # Calculate Dynamic Information
  if (length(object@Profile@adj_p_values) > 0) {
    sig.level <- object@Parameters@p_value
    nSigs <- length(object@Profile@adj_p_values[object@Profile@adj_p_values <= sig.level])
    if (all(object@Profile@adj_p_values > sig.level)) {
      cat("\nNo. of Significant Profiles: None found")
    } else {
      cat(paste("\nNo. of Significant Profiles:", nSigs, sep = " "))
    }
  }

  # Influential Genes if any
  if (ncol(object@Estimate@influential) > 0) {
    cat(paste("\nNo. of Influential Features:", ncol(object@Estimate@influential)))
  }
}

###############################################################################
#' @title Extract Information `.ScMaSigPro_show`
#'
#' @description
#' Extract details for the console cat of `.ScMaSigPro_show`.
#'
#' @param data A dataframe.
#' @param return_type Average Bin Size or Number of Bins.
#' @param bin_size_col Column name with the bins.
#' @param path_col Column name with the path information.
#'
#' @keywords internal
extract_info <- function(data, return_type = "avg_bin_size",
                         bin_size_col, path_col) {
  if (return_type == "avg_bin_size") {
    avg_sizes <- tapply(data[[bin_size_col]], data[[path_col]], mean)
    return(avg_sizes)
  } else if (return_type == "num_bins") {
    bin_counts <- table(data[[path_col]])
    return(bin_counts)
  } else {
    stop("Invalid return_type. Choose between 'avg_bin_size' and 'num_bins'.")
  }
}
