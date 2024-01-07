#' @title Get Features from scMaSigPro Object
#'
#' @description
#' This function extracts features from an object of class `scMaSigPro`
#' based on a specified query. The function query the coefficients matrix.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param rsq Coefficient of determination or R-squared value threshold.
#' @param includeInflu Whether to include genes with influential observations.
#' @param p_value Overall model significance.
#' @param query A string specifying the type of query. Valid options are 'unique' or 'union'.
#' Default is 'unique'.
#' @param change The effective change at the end of the pseudotime.
#' (Default is NULL)
#' @param strictly If `change != NULL`, selects features if all coefficients
#' have similar change.
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return A subset of matrix of coefficients.
#'
#' @seealso `scMaSigPro::sc.filter()`
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}.
#'
#' @export
queryCoeff <- function(scmpObj,
                       rsq = 0.7,
                       p_value = scmpObj@Parameters@p_value,
                       includeInflu = TRUE,
                       query = "pseudotime_path",
                       change = NULL,
                       strictly = FALSE,
                       verbose = TRUE) {
  # Check Validity of the object
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'ScMaSigPro'"
  )
  # Check for group_vector
  assert_that(!isEmpty(scmpObj@Design@groups.vector),
    msg = "'scmpObj@Design@groups.vector' is empty"
  )

  # Check for requested change
  if (!is.null(change)) {
    assert_that(
      all(
        change %in% c("increasing", "decreasing")
      ),
      msg = "Invalid change, please select on of 'increasing' or 'decreasing'"
    )
  }

  # Check query
  assert_that(
    all(
      query %in% c("pseudotime", "pseudotime_path", "path", "path_pseudotime")
    ),
    msg = "Invalid query, please select on of 'pseudotime', 'pseudotime_path' or 'path'"
  )

  # Check coefficents
  coeff_matrix <- showCoeff(scmpObj = scmpObj)

  # Apply R2 filter
  sol_df <- showSol(scmpObj = scmpObj)[, c(1, 2)] %>% as.data.frame()

  # Subset the sol
  sol_df <- sol_df[sol_df[["R-squared"]] >= rsq, , drop = FALSE]

  # Apply p-value filter
  sig_genes <- rownames(sol_df[sol_df[["p-value"]] <= p_value, , drop = FALSE])

  # Subset the coeff_matrix
  coeff_matrix <- coeff_matrix[sig_genes, , drop = FALSE]

  # All coeffients
  coeff_vector <- colnames(coeff_matrix)

  # Print Model formula
  message(paste("Polynomial-GLM Formula:", showPoly(scmpObj)))

  # Get groups
  compare_groups_vector <- unique(scmpObj@Design@groups.vector)

  # Generate group name vector
  avail_groups_vector <- unique(unlist(str_split(compare_groups_vector, "vs")))

  # Verbose
  if (verbose) {
    message(paste0(
      "Number of available groups ", length(compare_groups_vector),
      ", i.e. ", paste(compare_groups_vector, collapse = ", ")
    ))
  }

  # Query Simplification
  if (query == "pseudotime") {
    # Get temporal variable
    time_comp <- scmpObj@Parameters@bin_ptime_col

    # Get all the columns with time variable
    time_comp_vector <- grep(
      pattern = time_comp,
      x = coeff_vector,
      ignore.case = FALSE,
      value = TRUE
    )

    # Exclude any of the group effects
    group_effects <- paste(avail_groups_vector, collapse = "|")

    # Get all the columns with time only variable
    time_comp_vector <- grep(
      pattern = group_effects,
      x = time_comp_vector,
      ignore.case = FALSE,
      value = TRUE, invert = TRUE
    )

    # Get the coefficients
    time_coeff_matrix <- coeff_matrix[, time_comp_vector, drop = FALSE]

    # Use the apply function to check each row
    rows_with_all_nonzero <- apply(time_coeff_matrix, 1, function(row) all(row != 0))

    # Subset the data to keep only the rows where all values are non-zero
    time_coeff_matrix <- time_coeff_matrix[rows_with_all_nonzero, , drop = FALSE]

    # Check for requested change
    if (!is.null(change)) {
      if (change == "increasing") {
        if (strictly) {
          rows_increasing <- apply(time_coeff_matrix, 1, function(row) all(row > 0))
          time_coeff_matrix <- time_coeff_matrix[rows_increasing, , drop = FALSE]
        } else {
          time_coeff_matrix <- time_coeff_matrix[rowSums(time_coeff_matrix) > 0, , drop = FALSE]
        }
      } else if (change == "decreasing") {
        if (strictly) {
          rows_decreasing <- apply(time_coeff_matrix, 1, function(row) all(row < 0))
          time_coeff_matrix <- time_coeff_matrix[rows_decreasing, , drop = FALSE]
        } else {
          time_coeff_matrix <- time_coeff_matrix[rowSums(time_coeff_matrix) < 0, , drop = FALSE]
        }
      }
    }

    # Return
    return(time_coeff_matrix)
  } else if (query == "pseudotime_path" || query == "path_pseudotime") {
    # Get temporal variable
    time_comp <- scmpObj@Parameters@bin_ptime_col

    # Get all the columns with time variable
    time_comp_vector <- grep(
      pattern = time_comp,
      x = coeff_vector,
      ignore.case = FALSE,
      value = TRUE
    )

    # Exclude any of the group effects
    group_effects <- paste(avail_groups_vector, collapse = "|")

    # Get all the columns with time only variable
    time_comp_vector <- grep(
      pattern = group_effects,
      x = time_comp_vector,
      ignore.case = FALSE,
      value = TRUE, invert = FALSE
    )

    # Get the coeffients
    time_coeff_matrix <- coeff_matrix[, time_comp_vector, drop = FALSE]

    # Use the apply function to check each row
    rows_with_all_nonzero <- apply(time_coeff_matrix, 1, function(row) all(row != 0))

    # Subset the data to keep only the rows where all values are non-zero
    time_coeff_matrix <- time_coeff_matrix[rows_with_all_nonzero, , drop = FALSE]

    # If change is requested
    if (!is.null(change)) {
      if (change == "increasing") {
        if (strictly) {
          rows_increasing <- apply(time_coeff_matrix, 1, function(row) all(row > 0))
          time_coeff_matrix <- time_coeff_matrix[rows_increasing, , drop = FALSE]
        } else {
          time_coeff_matrix <- time_coeff_matrix[rowSums(time_coeff_matrix) > 0, , drop = FALSE]
        }
      } else if (change == "decreasing") {
        if (strictly) {
          rows_decreasing <- apply(time_coeff_matrix, 1, function(row) all(row < 0))
          time_coeff_matrix <- time_coeff_matrix[rows_decreasing, , drop = FALSE]
        } else {
          time_coeff_matrix <- time_coeff_matrix[rowSums(time_coeff_matrix) < 0, , drop = FALSE]
        }
      }
    }
    # Return
    return(time_coeff_matrix)
  }

  # Query Simplification
  else if (query == "path") {
    # Get temporal variable
    time_comp <- scmpObj@Parameters@bin_ptime_col

    # Exclude any of the group effects
    group_effects <- paste(avail_groups_vector, collapse = "|")

    # Get all the columns with time variable
    time_comp_vector <- grep(
      pattern = group_effects,
      x = coeff_vector,
      ignore.case = FALSE,
      value = TRUE
    )

    # Get all the columns with time only variable
    time_comp_vector <- grep(
      pattern = time_comp,
      x = time_comp_vector,
      ignore.case = FALSE,
      value = TRUE, invert = TRUE
    )

    # Get the coeffients
    time_coeff_matrix <- coeff_matrix[, time_comp_vector, drop = FALSE]

    # Use the apply function to check each row
    rows_with_all_nonzero <- apply(time_coeff_matrix, 1, function(row) all(row != 0))

    # Subset the data to keep only the rows where all values are non-zero
    time_coeff_matrix <- time_coeff_matrix[rows_with_all_nonzero, , drop = FALSE]

    # Check for requested change
    if (!is.null(change)) {
      if (change == "increasing") {
        if (strictly) {
          rows_increasing <- apply(time_coeff_matrix, 1, function(row) all(row > 0))
          time_coeff_matrix <- time_coeff_matrix[rows_increasing, , drop = FALSE]
        } else {
          time_coeff_matrix <- time_coeff_matrix[rowSums(time_coeff_matrix) > 0, , drop = FALSE]
        }
      } else if (change == "decreasing") {
        if (strictly) {
          rows_decreasing <- apply(time_coeff_matrix, 1, function(row) all(row < 0))
          time_coeff_matrix <- time_coeff_matrix[rows_decreasing, , drop = FALSE]
        } else {
          time_coeff_matrix <- time_coeff_matrix[rowSums(time_coeff_matrix) < 0, , drop = FALSE]
        }
      }
    }
    # Return
    return(time_coeff_matrix)
  }
}
