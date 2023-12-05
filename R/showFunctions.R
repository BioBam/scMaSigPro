## Show functions for scMaSigPro
## Description: Following functions are used to extract information from the scmpClass
## Object
# 1. showCoeff():
# 2. showInflu():
# 3. showTS():
# 4. showSol():
# 5. showSigProf():
# 6. .scmp_show(): For object cat
# 7. extract_info():
# 8. showParams():
# 9. showPoly():
# 10. showGroupCoeff():

###############################################################################

#' Show or Return the Coefficent matrix
#'
#' This function is used to view or return the coeffients of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param includeInflu description
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @export
showCoeff <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@sol) == c(0, 0)),
    msg = "Coeff is not computed yet"
  )

  # Extract
  coefficients <- scmpObj@scTFit@coefficients %>% as.data.frame()

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

#' Show or Return the matrix of influential genes
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @export
showInflu <- function(scmpObj, view = FALSE, return = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@influ.info) == c(0, 0)),
    msg = "tscore is not computed yet"
  )

  # Extract
  influ <- scmpObj@scTFit@influ.info

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

#' Show or Return the t scores
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param includeInflu logical, whether to add gene with influential data in the solution.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @export
showTS <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@t.score) == c(0, 0)),
    msg = "tscore is not computed yet"
  )

  # Extract
  tscore <- scmpObj@scTFit@t.score

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

#' Show or Return the Solution
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param includeInflu logical, whether to add gene with influential data in the solution.
#' @param return logical, whether to return the solution. If TRUE (default), returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @export
showSol <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@sol) == c(0, 0)),
    msg = "Sol is not computed yet"
  )

  # Extract
  sol <- scmpObj@scTFit@sol %>% as.data.frame()

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
#' Show or Return the Solution
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param includeInflu logical, whether to add gene with influential data in the solution.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#' @importFrom utils View
#' @export
showSigProf <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = FALSE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@sol) == c(0, 0)),
    msg = "Sol is not computed yet"
  )

  # Extract
  sol <- showSol(scmpObj, view = FALSE, return = TRUE, includeInflu = includeInflu) %>% as.data.frame()
  # Extract rownames
  bulk.counts <- scmpObj@compress.sce@assays@data@listData$bulk.counts
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
#' Show the terms of the polynomial term
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#'
#' @return Return the terms of the polynomial model.
#'
#' @export
showPoly <- function(scmpObj) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(all(!is.na(colnames(scmpObj@edesign@dis)) | length(colnames(scmpObj@edesign@dis) > 1)),
    msg = "Please setup the model first, using 'sc.make.design.matrix()'"
  )

  # Extract columns
  df.col <- colnames(scmpObj@edesign@dis)

  # Extract betas
  beta_names <- paste0("beta", seq(1:length(df.col)))

  # Generate formula string
  formula_parts <- vapply(seq_along(df.col),
    function(i) paste(beta_names[i], "*", df.col[i], sep = ""),
    FUN.VALUE = character(1)
  )

  # Make formula
  formula_string <- paste("beta0", paste(formula_parts, collapse = " + "), sep = " + ")

  return(formula_string)
}

###############################################################################
#' Show or Return the parameters used during the analysis
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#' @importFrom methods slot slotNames
#' @export
showParams <- function(scmpObj, view = FALSE, return = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Get slot names, assuming 'addParams' is a slot within 'scmpObj'
  all_slots <- slotNames(scmpObj)

  # Get 'addParams' slot data using the correct S4 accessor method
  addParamsData <- slot(scmpObj, "addParams")

  # Get all slots of 'addParams', assuming 'addParams' itself is an S4 object with slots
  params <- lapply(slotNames(addParamsData), function(parameter) {
    slot(addParamsData, parameter)
  })

  # Get the data
  params <- data.frame(
    parameters = slotNames(addParamsData),
    value = unlist(params)
  )

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
#' Show or Return the Group wise coefficents
#'
#' This function is used to view or return the group of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#' @param includeInflu description
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @export
showGroupCoeff <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@group.coeffs) == c(0, 0)),
    msg = "group.coeffs is not computed yet"
  )

  # Extract
  grpCoeff <- scmpObj@scTFit@group.coeffs %>% as.data.frame()

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
#' Show ScMaSigPro Object Information
#'
#' This method displays basic information about the ScMaSigPro object when the object
#' is printed in the console. The method is automatically called when the user writes
#' the name of the object in the console.
#'
#' @param object An object of class \code{scMaSigProClass}.
#'
#' @importFrom S4Vectors coolcat
#'
#' @keywords internal
#' @export
.scmp_show <- function(object) {
  # Show Basic information
  cat("Class: scmpClass\n")
  cat(paste0("nCells: ", ncol(object@sce), "\n"))
  cat(paste0("nFeatures: ", nrow(object@sce), "\n"))
  cat("Pseudotime Range:", paste(range(colData(object@sce)[[object@addParams@pseudotime_colname]])))

  # Calculate the Compression
  compressed.cell.metadata <- object@compress.sce@colData %>% as.data.frame()
  if (length(compressed.cell.metadata) > 0) {
    cat(paste("\nPaths:", paste(levels(as.factor(compressed.cell.metadata[[object@addParams@path_colname]])), collapse = ", ")))
    cat(paste0(
      "\nBinned Pseudotime: ", paste(range(compressed.cell.metadata[[object@addParams@bin_pseudotime_colname]]), collapse = "-"), "(Range), ",
      round(mean(compressed.cell.metadata[[object@addParams@bin_pseudotime_colname]]), 2), "(Mean), "
    ))

    # Extract info
    per_path_num_bin <- extract_info(compressed.cell.metadata, return_type = "num_bins", bin_size_col = object@addParams@bin_size_colname, object@addParams@path_colname)
    per_path_bin_size <- round(extract_info(compressed.cell.metadata, return_type = "avg_bin_size", bin_size_col = object@addParams@bin_size_colname, object@addParams@path_colname))

    # Paste
    cat("\nNumber of bins->", paste(names(per_path_num_bin), per_path_num_bin, sep = ": "))
    cat("\nAverage bin Size->", paste(names(per_path_bin_size), per_path_bin_size, sep = ": "))
  }

  # Influential Genes if any
  if (length(colnames(object@edesign@dis)) > 0) {
    cat(paste("\nPolynomial Order:", object@edesign@poly_degree))
  }

  # Calculate Dynamic Information
  if (length(object@scPVector@p.adjusted) > 0) {
    sig.level <- object@addParams@Q
    nSigs <- length(object@scPVector@p.adjusted[object@scPVector@p.adjusted <= sig.level])
    if (all(object@scPVector@p.adjusted > sig.level)) {
      cat("\nSig. Profiles (P-vector): None found")
    } else {
      cat(paste("\nSig. Models (sc.p.vector):", nSigs, sep = " "))
    }
  }

  # Influential Genes if any
  if (ncol(object@scTFit@influ.info) > 0) {
    cat(paste("\nNo. of Influential Features:", ncol(object@scTFit@influ.info)))
  }
}

# helper to extract the lineage info
extract_info <- function(data, return_type = "avg_bin_size", bin_size_col, path_col) {
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

#-X-#