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
#' @rdname showSol
#' @export
showCoeff <- function(scmpObj, includeInflu = TRUE) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check if the sol exist
    .check(!all(dim(scmpObj@Estimate@coefficient_matrix) == c(0, 0)),
        msg = "'coefficient_matrix' is not computed yet"
    )

    # Extract
    coefficients <- as.data.frame(scmpObj@Estimate@coefficient_matrix)

    if (!includeInflu) {
        influ.gene <- colnames(showInflu(scmpObj))
        coefficients <- coefficients[!(rownames(coefficients) %in% influ.gene), ]
    }

    return(coefficients)
}

###############################################################################
#' @rdname showSol
#' @export
showInflu <- function(scmpObj) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # After sc.t.fit(), a fit without influential genes has a 0 x 0 matrix
    .check(!all(dim(scmpObj@Estimate@significance_matrix) == c(0, 0)),
        msg = "'influential' is not computed yet"
    )

    # Extract
    influ <- scmpObj@Estimate@influential

    return(influ)
}

###############################################################################
#' @rdname showSol
#' @export
showTS <- function(scmpObj, includeInflu = TRUE) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check if the sol exist
    .check(!all(dim(scmpObj@Estimate@t_score_matrix) == c(0, 0)),
        msg = "tscore is not computed yet"
    )

    # Extract
    tscore <- scmpObj@Estimate@t_score_matrix

    if (!includeInflu) {
        influ.gene <- colnames(showInflu(scmpObj))
        tscore <- tscore[!(rownames(tscore) %in% influ.gene), ]
    }

    return(tscore)
}

###############################################################################
#' Results of the term selection
#'
#' Read the per-gene results that [sc.t.fit()] stores in a
#' [ScMaSigPro-class] object:
#' * `showSol()`: for each gene whose selected model has a significant
#'   coefficient, the p-value and R-squared of that model and the p-value of
#'   each coefficient of the full model (`NA` for the terms it leaves out).
#' * `showCoeff()`: the coefficients of these genes, 0 for the terms the
#'   model leaves out and `NA` for the terms that cannot be estimated.
#' * `showTS()`: their t statistics.
#' * `showGroupCoeff()`: the coefficients of the polynomial of each path.
#' * `showInflu()`: the genes with influential observations, one column per
#'   gene and one row per bin; an entry names the influence measures that
#'   flag the bin (DFFITS and Cook's distance, see
#'   [stats::influence.measures()]).
#' * `showSigProf()`: the pseudo-bulk counts of the genes of `showSol()`.
#'
#' @param includeInflu Whether to keep the genes with influential
#'   observations, the columns of `showInflu()`.
#'
#' @return `showSol()`, `showCoeff()` and `showGroupCoeff()`: a data frame
#' with one row per gene. `showTS()`: a numeric matrix with one row per gene.
#' `showInflu()`: a character matrix with one row per bin and one column per
#' gene, 0 x 0 when no gene has influential observations. `showSigProf()`: a
#' `dgCMatrix` with one row per gene and one column per bin.
#'
#' @family results
#'
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' head(showSol(scmp.ob))
#' head(showGroupCoeff(scmp.ob))
#' colnames(showInflu(scmp.ob))
#'
#' @export
showSol <- function(scmpObj, includeInflu = TRUE) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check if the sol exist
    .check(!all(dim(scmpObj@Estimate@significance_matrix) == c(0, 0)),
        msg = "'significance_matrix' is not computed yet"
    )

    # Extract
    sol <- as.data.frame(scmpObj@Estimate@significance_matrix)

    if (!includeInflu) {
        influ.gene <- colnames(showInflu(scmpObj))
        sol <- sol[!(rownames(sol) %in% influ.gene), ]
    }

    return(sol)
}

###############################################################################
#' @rdname showSol
#' @export
showSigProf <- function(scmpObj, includeInflu = FALSE) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check if the sol exist
    .check(!all(dim(scmpObj@Estimate@significance_matrix) == c(0, 0)),
        msg = "Sol is not computed yet"
    )

    # Extract
    sol <- as.data.frame(showSol(scmpObj, includeInflu = includeInflu))
    # Extract rownames
    .check_dense_counts(scmpObj)
    bulk.counts <- assay(scmpObj@Dense, "bulk.counts", withDimnames = FALSE)
    bulk.counts <- bulk.counts[rownames(bulk.counts) %in% rownames(sol), , drop = FALSE]

    if (!includeInflu) {
        influ.gene <- colnames(showInflu(scmpObj))
        bulk.counts <- bulk.counts[!(rownames(bulk.counts) %in% influ.gene), ]
    }

    return(bulk.counts)
}

###############################################################################
#' @rdname showParams
#' @export
showPoly <- function(scmpObj) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check if the sol exist
    .check(
        all(!is.na(colnames(scmpObj@Design@predictor_matrix)) | length(
            colnames(scmpObj@Design@predictor_matrix) > 1
        )),
        msg = "Please setup the model first, using 'sc.set.poly()'"
    )

    # Extract columns
    df.col <- colnames(scmpObj@Design@predictor_matrix)

    # Extract betas
    beta_names <- paste0("beta", seq_along(df.col))

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
#' Parameters and model formula of an analysis
#'
#' `showParams()` lists the parameters that the steps of the workflow store
#' in a [ScMaSigPro-class] object (see [ParameterConfig-class]).
#' `showPoly()` writes the full model that [sc.set.poly()] set up, with one
#' coefficient `beta1`, `beta2`, ... per column of the design matrix.
#'
#' @return `showParams()`: a data frame with the columns `parameters` and
#' `value`, one row per parameter. `showPoly()`: a character string.
#'
#' @family results
#'
#' @importFrom methods slot slotNames
#'
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' head(showParams(scmp.ob))
#' showPoly(scmp.ob)
#'
#' @export
showParams <- function(scmpObj) {
    # Check Object Validity
    .check_scmp(scmpObj)

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
    distribution <- scmpObj@Parameters@distribution
    distribution.f <- data.frame(
        parameters = "distribution",
        value = if (identical(distribution$theta_estimate, "gene")) {
            "Negative Binomial(theta per gene)"
        } else {
            paste0(distribution[["family"]])
        }
    )


    # Add
    params <- rbind(
        params.frame,
        distribution.f
    )
    rownames(params) <- NULL

    return(params)
}

###############################################################################
#' @rdname showSol
#' @export
showGroupCoeff <- function(scmpObj, includeInflu = TRUE) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # Check if the sol exist
    .check(!all(dim(scmpObj@Estimate@path_coefficient_matrix) == c(0, 0)),
        msg = "path_coefficient_matrix is not computed yet"
    )

    # Extract
    grpCoeff <- as.data.frame(scmpObj@Estimate@path_coefficient_matrix)

    if (!includeInflu) {
        influ.gene <- colnames(showInflu(scmpObj))
        grpCoeff <- grpCoeff[!(rownames(grpCoeff) %in% influ.gene), ]
    }

    return(grpCoeff)
}

###############################################################################
#' @title Print a ScMaSigPro object
#'
#' @description
#' `.ScMaSigPro_show()` prints a summary of a [ScMaSigPro-class] object; the
#' `show()` method calls it.
#'
#' @param object A [ScMaSigPro-class] object.
#'
#' @keywords internal
#'
#' @seealso [ScMaSigPro-class]
#'
#' @return `NULL`, invisibly; called for what it prints.
.ScMaSigPro_show <- function(object) {
    # Show Basic information
    cat("Class: ScMaSigProClass\n")
    cat(paste0("nCells: ", ncol(object@Sparse), "\n"))
    cat(paste0("nFeatures: ", nrow(object@Sparse), "\n"))
    cat("Pseudotime Range:", paste(round(
        range(SingleCellExperiment::colData(object@Sparse)[[object@Parameters@ptime_col]]), 3
    )))
    # Once binned, the paths of the bins (sc.restruct() folds the root and
    # link paths into the end paths), in the order of the cells
    path_col <- object@Parameters@path_col
    paths <- unique(SingleCellExperiment::colData(object@Sparse)[[path_col]])
    bin_paths <- unique(cDense(object)[[path_col]])
    if (length(bin_paths) > 0) {
        paths <- c(intersect(paths, bin_paths), setdiff(bin_paths, paths))
    }
    cat(paste("\nBranching Paths:", paste(paths, collapse = ", ")))

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

    # Saved before the version slot or the count filter parameters existed,
    # or with the log_offset parameter
    if (!.scmp_current(object)) {
        cat(paste(
            "\nSaved by an earlier version of scMaSigPro: update it with",
            "updateObject()."
        ))
    }
}

###############################################################################
#' @title Bins per path for the printed summary
#'
#' @description
#' `extract_info()` gives the mean bin size or the number of bins of each path
#' for `.ScMaSigPro_show()`.
#'
#' @param data The bin metadata.
#' @param return_type `"avg_bin_size"` or `"num_bins"`.
#' @param bin_size_col The column with the number of cells of each bin.
#' @param path_col The column with the path labels.
#'
#' @keywords internal
#'
#' @seealso [ScMaSigPro-class], whose `show()` method prints it.
#'
#' @return One value per path: the mean bin size (`"avg_bin_size"`) or the
#' number of bins (`"num_bins"`, a table).
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
