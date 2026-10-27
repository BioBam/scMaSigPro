#' @title Query the coefficients of the selected models
#'
#' @description
#' `queryCoeff()` returns the coefficients of a group of terms for the genes
#' whose selected model ([sc.t.fit()]) passes `rsq` and `p_value` and has a
#' non-zero coefficient for every term of the group. It also prints the
#' formula of the model.
#'
#' @param query The terms: `"pseudotime"`, the powers of pseudotime of the
#'   reference path; `"pseudotime_path"` or `"path_pseudotime"`, their products
#'   with the path dummies; `"path"`, the path dummies.
#' @param change `NULL`, or `"increasing"` or `"decreasing"`: keep the genes
#'   whose coefficients of the query add up to a positive (negative) number,
#'   or with `strictly = TRUE` are all positive (negative).
#' @param strictly See `change`.
#'
#' @return A data frame of coefficients: one row per gene that passes the
#' thresholds and the query, one column per queried term.
#'
#' @seealso [showCoeff()] for the coefficients of every gene.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}.
#'
#' @family results
#' @inheritParams scmp-params
#' @inheritParams sc.filter
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' # Genes whose coefficients of pseudotime are all positive
#' queryCoeff(scmp.ob,
#'     query = "pseudotime", change = "increasing", strictly = TRUE,
#'     verbose = FALSE
#' )
#'
#' @export
queryCoeff <- function(scmpObj,
                       rsq = 0.7,
                       p_value = scmpObj@Parameters@p_value,
                       query = "pseudotime_path",
                       change = NULL,
                       strictly = FALSE,
                       verbose = TRUE) {
    # Check Validity of the object
    .check_scmp(scmpObj)
    # Check for group_vector
    .check(!isEmpty(scmpObj@Design@groups.vector),
        msg = "'scmpObj@Design@groups.vector' is empty"
    )

    # Check for requested change
    if (!is.null(change)) {
        change <- .check_choice(change, c("increasing", "decreasing"),
            msg = paste(
                "Invalid change, please select one of 'increasing' or",
                "'decreasing'"
            )
        )
    }

    # Check query
    query <- .check_choice(query,
        c("pseudotime", "pseudotime_path", "path", "path_pseudotime"),
        msg = paste(
            "Invalid query, please select one of 'pseudotime',",
            "'pseudotime_path' or 'path'"
        )
    )

    # Check coefficents
    coeff_matrix <- showCoeff(scmpObj = scmpObj)

    # Apply R2 filter
    sol_df <- as.data.frame(showSol(scmpObj = scmpObj)[, c(1, 2)])

    # Subset the sol
    sol_df <- sol_df[sol_df[["R-squared"]] >= rsq, , drop = FALSE]

    # Apply p-value filter
    sig_genes <- rownames(sol_df[sol_df[["p-value"]] <= p_value, , drop = FALSE])

    # Subset the coeff_matrix
    coeff_matrix <- coeff_matrix[sig_genes, , drop = FALSE]

    # All coeffients
    coeff_vector <- colnames(coeff_matrix)

    # Print Model formula
    message("Polynomial-GLM Formula: ", showPoly(scmpObj))

    # Get groups
    compare_groups_vector <- unique(scmpObj@Design@groups.vector)

    # Generate group name vector
    avail_groups_vector <- unique(unlist(
        strsplit(compare_groups_vector, "vs", fixed = TRUE)
    ))

    # Verbose
    if (verbose) {
        message(
            "Number of available groups ", length(compare_groups_vector),
            ", i.e. ", paste(compare_groups_vector, collapse = ", ")
        )
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
