#################################################
## Author: Priyansh Srivastava ##################
## Year: 2023 ###################################
## Email: spriyansh29@gmail.com #################
#################################################
# Internal function for pseudo-bulking counts along Pseudotime. They are called
# by sc.squeeze()
###############################################################################
#' @title Number of cells of a bin
#'
#' @description
#' `calc_bin_size()` counts the cells of a bin: the names separated by `|`
#' in its `clus_mem_col` column.
#'
#' @param x A row of the bin metadata.
#' @param clus_mem_col The name of the column with the cells of the bin.
#'
#' @return The number of cells of the bin.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @keywords internal
#' @family binning helpers
calc_bin_size <- function(x, clus_mem_col = "scmp_cluster_members") {
    # Split the bin's cell members at each "|" (a fixed string, not a regular
    # expression); the number of pieces is the number of cells in the bin
    size <- length(strsplit(x[[clus_mem_col]], "|", fixed = TRUE)[[1]])

    # Convert the 'size' variable to a numeric value and return it as the result of the function
    return(as.numeric(size))
}
###############################################################################
#' @title Number of bins of a histogram rule
#'
#' @description
#' `estBinSize()` gives the number of bins of one of the histogram rules
#' "Sturges", "Sqrt", "Rice", "Doane", "Freedman.Diaconis" and "Scott.Normal"
#' for a pseudotime vector, multiplied by `drop_fac`.
#'
#' @param time_vector The pseudotime values to bin.
#' @param nPoints The number of values of `time_vector`.
#' @param drop_fac The factor that multiplies the number of bins of the rule;
#'   values below 1 give fewer, wider bins.
#' @param bin_method The rule (see Details).
#'
#' @details
#' With n the number of points (`nPoints`), the rules give:
#' \describe{
#'   \item{"Sturges"}{log2(n) + 1 bins (Sturges 1926). Srivastava et al.
#'   (2024), Supplementary Section S1.1, writes the number of bins as
#'   \eqn{B = [(\log_2(N)+1) \times k]}{B = [(log2(N)+1) * k]}, with k the
#'   drop factor of their Section 2.1, and every other rule the same way,
#'   with the brackets around the whole product: the rule times `drop_fac`
#'   is rounded up.}
#'   \item{"Sqrt"}{sqrt(n) bins, a common rule of thumb.}
#'   \item{"Rice"}{2 n^(1/3) bins (Lane and Scott), close to the
#'   oversmoothed number of bins (2n)^(1/3) of Terrell and Scott (1985).}
#'   \item{"Doane"}{1 + log2(n) + log2(1 + |g| / s) bins, with g the sample
#'   skewness of the points and s = sqrt(6 (n - 2) / ((n + 1) (n + 3)))
#'   (Doane 1976).}
#'   \item{"Freedman.Diaconis"}{bins of width 2 IQR / n^(1/3), with IQR the
#'   interquartile range of the points (Freedman and Diaconis 1981).}
#'   \item{"Scott.Normal"}{bins of width 3.49 sd / n^(1/3), with sd the
#'   standard deviation of the points (Scott 1979).}
#' }
#' The widths of Freedman-Diaconis and Scott are turned into a number of bins
#' by dividing the range of the points by them. A width of 0 (all points
#' equal, or for Freedman-Diaconis equal first and third quartiles) or an
#' undefined width (Scott with a single point) gives one bin, with a message.
#'
#' @return The number of bins of the rule times `drop_fac`, not rounded;
#' [sc.squeeze()] rounds it up.
#'
#' @references
#' Sturges HA (1926). The choice of a class interval. Journal of the
#' American Statistical Association 21(153), 65-66.
#' \doi{10.1080/01621459.1926.10502161}
#'
#' Doane DP (1976). Aesthetic frequency classifications. The American
#' Statistician 30(4), 181-183. \doi{10.1080/00031305.1976.10479172}
#'
#' Scott DW (1979). On optimal and data-based histograms. Biometrika 66(3),
#' 605-610. \doi{10.1093/biomet/66.3.605}
#'
#' Freedman D, Diaconis P (1981). On the histogram as a density estimator:
#' L2 theory. Zeitschrift für Wahrscheinlichkeitstheorie und verwandte
#' Gebiete 57(4), 453-476. \doi{10.1007/BF01025868}
#'
#' Terrell GR, Scott DW (1985). Oversmoothed nonparametric density
#' estimates. Journal of the American Statistical Association 80(389),
#' 209-214. \doi{10.1080/01621459.1985.10477163}
#'
#' Lane DM, Scott D. Histograms. In Online Statistics Education: A
#' Multimedia Course of Study, Rice University.
#' <https://onlinestatbook.com/2/graphing_distributions/histograms.html>
#'
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom stats IQR sd
#'
#' @keywords internal
#' @family binning helpers
estBinSize <- function(time_vector, nPoints, drop_fac, bin_method) {
    estBins <- switch(bin_method,
        "Freedman.Diaconis" = {
            # Freedman-Diaconis rule: bins of width twice the interquartile
            # range over the cube root of the number of points
            width <- 2 * IQR(time_vector) / nPoints^(1 / 3)
            .width_to_count(time_vector, width, bin_method)
        },
        "Sqrt" = {
            # Square root rule: the square root of the number of points
            nPoints^(1 / 2)
        },
        "Sturges" = {
            # Sturges' rule: the log (base 2) of the number of points, plus 1
            log2(nPoints) + 1
        },
        "Rice" = {
            # Rice rule: twice the cube root of the number of points
            2 * nPoints^(1 / 3)
        },
        "Doane" = {
            # Doane's formula: Sturges' rule corrected for the skewness
            sigma <- ((6 * (nPoints - 2)) / ((nPoints + 1) * (nPoints + 3)))^(1 / 2)
            sk <- .skewness(time_vector)
            1 + log2(nPoints) + log2(1 + (abs(sk) / sigma))
        },
        "Scott.Normal" = {
            # Scott's normal reference rule: bins of width 3.49 standard
            # deviations over the cube root of the number of points
            width <- 3.49 * sd(time_vector) / nPoints^(1 / 3)
            .width_to_count(time_vector, width, bin_method)
        },
        stop(
            "Invalid bin_method: ", bin_method, ". Please choose one of the ",
            "following: 'Freedman.Diaconis', 'Sqrt', 'Sturges', 'Rice', ",
            "'Doane', 'Scott.Normal'"
        )
    )


    # Scale the estimated number of bins by the drop factor.
    estBins <- drop_fac * estBins

    return(estBins)
}
###############################################################################
#' @title Number of Bins of a Given Width
#'
#' @description
#' `.width_to_count()` turns the bin width of the Freedman-Diaconis or the
#' Scott rule into a number of bins: the range of `x` divided by the width.
#' A width of 0 (an interquartile range or a standard deviation of 0) or an
#' undefined width (a single point) gives one bin, with a message.
#'
#' @param x A numeric vector: the points to bin.
#' @param width The bin width given by the rule.
#' @param bin_method The rule, "Freedman.Diaconis" or "Scott.Normal".
#'
#' @return A number, not rounded; `NA` when `x` has missing or infinite
#' values, which `sc.squeeze()` cannot bin with any rule.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @noRd
.width_to_count <- function(x, width, bin_method) {
    if (!all(is.finite(x))) {
        return(NA_real_)
    }
    if (isTRUE(width > 0)) {
        return(diff(range(x)) / width)
    }
    spread <- if (bin_method == "Freedman.Diaconis") {
        "interquartile range"
    } else {
        "standard deviation"
    }
    message(
        "The ", spread, " of the ", length(x), " pseudotime value(s) is ",
        if (is.na(width)) "undefined" else "0", ", so bin_method = '",
        bin_method, "' gives no bin width; they are counted as one bin."
    )
    1
}
###############################################################################
#' @title Sample Skewness
#'
#' @description
#' `.skewness()` is the sample skewness b1 = m3 / s^3 of `x`, where m3 is the
#' third central moment and s the standard deviation: type 3 of Joanes and
#' Gill (1998), the default of `e1071::skewness()`, which it replaces. It is
#' computed as sqrt(n) * sum(d^3) / sum(d^2)^(3/2) * (1 - 1/n)^(3/2) with
#' d = x - mean(x), the same operations in the same order, so the Doane rule
#' of `estBinSize()` gets the same number.
#'
#' @param x A numeric vector.
#'
#' @return A number: `NA` if `x` has missing values, `NaN` if all its values
#' are equal.
#'
#' @references Joanes, D. N. and Gill, C. A. (1998). Comparing measures of
#' sample skewness and kurtosis. The Statistician, 47, 183-189.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @noRd
.skewness <- function(x) {
    if (anyNA(x)) {
        return(NA)
    }
    n <- length(x)
    d <- x - mean(x)
    sqrt(n) * sum(d^3) / sum(d^2)^(3 / 2) * (1 - 1 / n)^(3 / 2)
}
###############################################################################
#' @title Summaries of the selected model of a gene
#'
#' @description
#' `extract_fitting()` collects, for one gene, the p-value and R-squared of
#' its selected model, the coefficients, p-values and t statistics of its
#' terms and its influential observations; [sc.t.fit()] calls it for each
#' gene.
#'
#' @param reg The selected model, a `glm` fit.
#' @param lmf The full model.
#' @param model.glm.0 The intercept-only model.
#' @param dis The design, a data frame with one column per term.
#' @param family The family of the fits.
#' @param name The name of the gene: the row name of `sol` and the column
#'   name of its entry in `influ.info`.
#' @param vars.in The names of the terms (the columns of `dis`); they fix the
#'   order of the coefficients, p-values and t statistics.
#' @param alfa The level at which a coefficient is significant.
#' @param influ.info The matrix of influential observations so far; a column
#'   for this gene is added when it has influential observations.
#'
#' @keywords internal
#'
#' @seealso [sc.t.fit()], which calls it.
#'
#' @return A list: the p-value (`p_value`) and R-squared (`bondad`) of the
#' selected model; the p-values (`p_valor`), coefficients (`coeff`) and t
#' statistics (`t`) of every term; the counts of the gene (`sig_profiles`);
#' its row of the results (`sol`, `NULL` when no coefficient is
#' significant); the updated `influ.info`; and the gene name
#' (`feature_name`).
extract_fitting <- function(reg, lmf, model.glm.0, dis, family, name, vars.in, alfa, influ.info) {
    sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
    y <- reg$y
    result <- summary(lmf)
    novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 4]))]
    influ <- influence.measures(reg)$is.inf
    influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
    influ1 <- which(apply(influ, 1, all))
    if (length(influ1) != 0) {
        paste.names <- function(a) {
            paste(names(a)[a], collapse = "/")
        }
        match <- match(rownames(dis), rownames(influ))
        influ <- as.data.frame(apply(influ, 1, paste.names))
        influ.info <- cbind(influ.info, influ[match, ])
        colnames(influ.info)[ncol(influ.info)] <- name
        influ.info <- as.matrix(influ.info)
    }
    result <- summary(reg)
    gaussian_ok <- !(result$aic == -Inf) && !is.na(result$aic) &&
        family$family == "gaussian"
    if (gaussian_ok || family$family != "gaussian") {
        if (family$family == "gaussian") {
            test <- anova(model.glm.0, reg, test = "F")
            p.value <- test[6][2, 1]
        } else {
            test <- anova(model.glm.0, reg, test = "Chisq")
            p.value <- test[5][2, 1]
        }
        # Computing goodness of fitting:
        # Why: the R-squared of a GLM is the share of the null deviance
        # that the model explains; a model worse than the intercept counts
        # as 0.
        bondad <- (reg$null.deviance - reg$deviance) / reg$null.deviance
        if (bondad < 0) {
            bondad <- 0
        }
        beta.coeff <- result$coefficients[, 1]
        beta.p.valor <- result$coefficients[, 4]
        coeff <- rep(0, (length(vars.in) + 1))
        # Terms that are not estimable in the full model get NA: every column
        # of dis with a name in novar (a name that is absent matches none)
        coeff[which(colnames(dis) %in% novar) + 1] <- NA
        p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 1)))

        if (result$coefficients[, 4][
            rownames(result$coefficients) == "(Intercept)"
        ] < alfa) {
            coeff[1] <- result$coefficients[, 1][
                rownames(result$coefficients) == "(Intercept)"
            ]
            p.valor[1] <- result$coefficients[, 4][
                rownames(result$coefficients) == "(Intercept)"
            ]
            t[1] <- result$coefficients[, 3][
                rownames(result$coefficients) == "(Intercept)"
            ]
        }
        for (j in 2:length(coeff)) {
            if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
                coeff[j] <- result$coefficients[, 1][
                    rownames(result$coefficients) == vars.in[j - 1]
                ]
                p.valor[j] <- result$coefficients[, 4][
                    rownames(result$coefficients) == vars.in[j - 1]
                ]
                t[j] <- result$coefficients[, 3][
                    rownames(result$coefficients) == vars.in[j - 1]
                ]
            }
        }
        if (!all(is.na(p.valor))) {
            sol <- rbind(sol, as.numeric(c(
                p.value, bondad,
                p.valor
            )))
            coefficients <- rbind(coefficients, coeff)
            t.score <- rbind(t.score, t)
            sig.profiles <- rbind(sig.profiles, y)
            h <- nrow(sol)
            rownames(sol)[h] <- name
            rownames(coefficients)[h] <- name
            rownames(t.score)[h] <- name
            rownames(sig.profiles)[h] <- name
        }
    }
    # Return Calculation
    return(list(
        p_value = p.value,
        bondad = bondad,
        p_valor = p.valor,
        coeff = coeff,
        t = t,
        sig_profiles = y,
        sol = sol,
        influ.info = influ.info,
        feature_name = name
    ))
}

###############################################################################
#' @title Bins of equal width
#'
#' @description
#' `extract_interval()` divides the range of `time.vector` into `nBins` bins
#' of equal width, whose bounds are the breaks
#' `seq(min, max, length.out = nBins + 1)`, and counts the points of each bin.
#'
#' @param time.vector The pseudotime values.
#' @param nBins The number of bins.
#' @param bin.size The name of the column for the number of points.
#' @param lbound,ubound The names of the columns for the lower and upper
#'   bounds.
#'
#' @return A data frame with one row per bin, in increasing order, and the
#' columns `lbound`, `ubound` and `bin.size`.
#'
#' @keywords internal
#' @family binning helpers
extract_interval <- function(time.vector, nBins = 1, bin.size, lbound,
                             ubound) {
    # Count the points in nBins bins of equal width spanning the range, with
    # the intervals of .assign_bins(). The bounds are the numeric breaks, not
    # the labels of cut(), which round them to 3 significant digits
    time_range <- range(time.vector)
    breaks <- seq(time_range[1], time_range[2], length.out = nBins + 1)
    codes <- cut(time.vector,
        breaks = breaks, right = FALSE, include.lowest = TRUE, labels = FALSE
    )
    new_bin_table_current <- data.frame(
        breaks[-length(breaks)], breaks[-1],
        as.numeric(tabulate(codes, nbins = length(breaks) - 1L))
    )
    colnames(new_bin_table_current) <- c(lbound, ubound, bin.size)
    return(new_bin_table_current)
}

###############################################################################
#' @title Assign Points to Bins
#'
#' @description
#' `.assign_bins()` gives each point the bin whose interval `[lower, upper)`
#' holds it, the last bin being closed on both ends, as in Srivastava et al.
#' (2024), Section 2.1: a point on the bound of two bins goes to the upper
#' one. On contiguous bins, these are the intervals of
#' `cut(right = FALSE, include.lowest = TRUE)`.
#'
#' @param x Numeric vector of points (pseudotime values).
#' @param lower,upper Numeric vectors: the bounds of the bins, in increasing
#' order.
#'
#' @return An integer vector: the bin of each point (its position in `lower`),
#' `NA` for a point that no bin holds.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @noRd
.assign_bins <- function(x, lower, upper) {
    idx <- rep(NA_integer_, length(x))
    last <- length(lower)
    for (j in seq_along(lower)) {
        below <- if (j == last) x <= upper[j] else x < upper[j]
        idx[x >= lower[j] & below] <- j
    }
    idx
}

###############################################################################
#' @title Split the bins that are too large
#'
#' @description
#' `optimize_bin_max()` splits every bin that holds more than `max_allowed`
#' points into the number of bins that `method` and `drop` give for its
#' points ([estBinSize()]), of equal width over their range; the first new
#' bin starts at the lower bound of the bin and the last ends at its upper
#' bound, so the new bins cover it without gaps. The other bins are kept.
#'
#' @param bin_table A data frame of bins, one row per bin, with the columns
#'   `lbound`, `ubound` and `bin.size`.
#' @param max_allowed The largest number of points of a bin that is kept.
#' @param verbose Whether to print a message for each bin.
#' @param time_vector The pseudotime of the cells of the path.
#' @param lbound,ubound The names of the columns of the lower and upper
#'   bounds.
#' @param bin.size The name of the column of the number of points.
#' @param method The rule that gives the number of new bins.
#' @param drop The factor that multiplies that number (`drop_fac`).
#'
#' @details
#' [sc.squeeze()] passes as `max_allowed` the maximum of Srivastava et al.
#' (2024), Section 2.2: the mean plus the standard deviation of the sizes of
#' the bins of the path, set once from the bins of the rule. "Bins larger
#' than the maximum size are divided", and "The process is repeated until all
#' bins conform to the established size criteria": [sc.squeeze()] calls
#' `optimize_bin_max()` again, with the same `max_allowed`, until no bin is
#' larger.
#'
#' @return A data frame with the columns of `bin_table`, one row per bin,
#' where each split bin is replaced by its parts.
#' @keywords internal
#' @family binning helpers
optimize_bin_max <- function(bin_table, max_allowed, verbose = TRUE,
                             time_vector, lbound, ubound, bin.size, method,
                             drop) {
    # Initiate an empty uniform bin
    uniform_bin_df <- data.frame(matrix(NA, ncol = ncol(bin_table)))
    colnames(uniform_bin_df) <- colnames(bin_table)

    # The bin of every point, by the intervals of .assign_bins()
    bin_of <- .assign_bins(
        time_vector, bin_table[[lbound]], bin_table[[ubound]]
    )

    # Binning one-by one
    for (i in seq_len(nrow(bin_table))) {
        # histo_bin
        current_bin_df <- c(unlist(bin_table[i, , drop = FALSE]))

        # Get size of the bins
        current_bin_size <- current_bin_df[[bin.size]]

        # Check if the bin size is big
        if (current_bin_size > max_allowed) {
            # get pseudotime
            pTime.for.big.interval <- time_vector[which(bin_of == i)]

            # Estimate how many splits are required
            potential_splits <- ceiling(estBinSize(
                time_vector = pTime.for.big.interval, nPoints = length(pTime.for.big.interval),
                drop_fac = drop, bin_method = method
            ))

            # Run the extraction
            small_splitted_bin_df <- extract_interval(
                time.vector = pTime.for.big.interval,
                nBins = potential_splits,
                bin.size = bin.size, lbound = lbound, ubound = ubound
            )
            # The new bins span the range of the cells; the outer two take
            # the bounds of the bin, so that the bins of the path follow
            # each other without gaps
            n_new <- nrow(small_splitted_bin_df)
            small_splitted_bin_df[[lbound]][1] <- current_bin_df[[lbound]]
            small_splitted_bin_df[[ubound]][n_new] <- current_bin_df[[ubound]]
            # Validation
            if (verbose) {
                message(
                    "Splitting bin with ", current_bin_size, " cells into ",
                    nrow(small_splitted_bin_df), " with sizes as ",
                    paste(small_splitted_bin_df[[bin.size]], collapse = ", ")
                )
            }

            # Add new Row
            uniform_bin_df <- rbind(uniform_bin_df, small_splitted_bin_df)
        } else if (current_bin_size <= max_allowed) {
            uniform_bin_df <- rbind(uniform_bin_df, current_bin_df)
            if (verbose) {
                if (current_bin_size != 0) {
                    message("Skipping bin with size ", current_bin_size)
                }
            }
        }
    }
    uniform_bin_df <- uniform_bin_df[-1, ]

    # Correct the rows
    rownames(uniform_bin_df) <- NULL

    # Return the new data
    return(uniform_bin_df)
}

################################################################################
#' @title Pseudo-bulk counts of the bins
#'
#' @description
#' `pb_counts()` adds up, or averages, the counts of the cells of each bin
#' into the pseudo-bulk counts (`bulk.counts`) of the `Dense` slot.
#' [sc.squeeze()] and [sc.restruct()] call it; call it after changing the
#' bins with [cDense()].
#'
#' @param cluster_count_by How the counts of the cells of a bin are combined:
#'   `"sum"` or `"mean"` (rounded to whole numbers).
#'
#' @details
#' The bins are the rows of the bin metadata ([cDense()]): their row names,
#' or the `bin_col` column when there are none, name the columns of the
#' counts, and `bin_mem_col` gives their cells.
#'
#' @return The object whose `Dense` slot holds the pseudo-bulk counts
#' (`bulk.counts`, genes by bins) and the bin metadata. The offsets of
#' the `Design` slot and the `Profile`, `Estimate` and
#' `Significant` slots, computed from earlier counts, are emptied; the
#' design is kept.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @family binning helpers
#' @inheritParams scmp-params
#'
#' @references
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' # The mean instead of the sum of the counts of the cells of each bin
#' scmp <- pb_counts(scmp.ob, cluster_count_by = "mean")
#' eDense(scmp)[1:3, 1:4]
#'
#' @export

pb_counts <- function(scmpObj,
                      bin_mem_col = scmpObj@Parameters@bin_mem_col,
                      bin_col = scmpObj@Parameters@bin_col,
                      assay_name = "counts",
                      cluster_count_by = "sum") {
    # Check Object Validity
    .check_scmp(scmpObj)
    .check_strings(
        bin_mem_col = bin_mem_col, bin_col = bin_col, assay_name = assay_name
    )
    cluster_count_by <- .check_choice(cluster_count_by, c("sum", "mean"),
        msg = "Invalid 'aggregate' value. Please choose either 'mean' or 'sum'."
    )

    # Count slot
    .check(
        all(
            assay_name %in% assayNames(scmpObj@Sparse)
        ),
        msg = paste0("'", assay_name, "' ", "doesn't exist in scmpObj.")
    )

    # Get assay
    counts <- assay(scmpObj@Sparse, assay_name, withDimnames = FALSE)

    # Get Pseudobulk Profile
    bin_data <- colData(scmpObj@Dense)
    pseudo_bulk_profile <- as.data.frame(bin_data)

    .check(bin_mem_col %in% colnames(pseudo_bulk_profile),
        msg = paste0("'", bin_mem_col, "' does not exist in level.meta.data")
    )
    .check(bin_col %in% colnames(pseudo_bulk_profile),
        msg = paste0("'", bin_col, "' does not exist in level.meta.data")
    )

    # Bins are named by the rownames of the bin metadata; these become the
    # column names of the assay (sc.set.poly() and sc.p.vector() rely on it)
    bin_ids <- rownames(bin_data)
    if (is.null(bin_ids)) {
        bin_ids <- as.character(bin_data[[bin_col]])
        rownames(bin_data) <- bin_ids
    }
    .check(anyDuplicated(bin_ids) == 0,
        msg = paste(
            "Bin names (the rownames of the 'Dense' cell metadata) must be",
            "unique."
        )
    )

    # Get the meta-information for pseudobulking
    meta.info <- pseudo_bulk_profile[, c(bin_mem_col, bin_col)]

    # Pseudo-bulk the counts of each bin
    pb.counts <- lapply(seq_len(nrow(meta.info)), function(i) {
        # Get the bin.info
        bin <- meta.info[i, , drop = FALSE]

        # Split the row
        cell.vector <- strsplit(bin[[1]], "|", fixed = TRUE)[[1]]

        # Get col cells
        col_indices <- which(colnames(counts) %in% cell.vector)

        # Subset the matrix using these indices
        bin_matrix <- as.matrix(counts[, col_indices, drop = FALSE])

        # Get Pseudobulked-counts
        pb.vector <- switch(cluster_count_by,
            "mean" = as.matrix(round(rowMeans(bin_matrix))),
            "sum"  = as.matrix(rowSums(bin_matrix)),
            stop("Invalid 'aggregate' value. Please choose either 'mean' or 'sum'.")
        )

        # Return
        return(pb.vector)
    })

    # Bind the per-bin columns into a matrix and set the row names
    pb.counts <- do.call(cbind, pb.counts)
    rownames(pb.counts) <- rownames(counts)
    colnames(pb.counts) <- bin_ids

    # Rebuild Dense so that the assay and the bin metadata are parallel
    scmpObj@Dense <- SingleCellExperiment(
        assays = list(bulk.counts = as(pb.counts, "dgCMatrix")),
        colData = bin_data
    )
    scmpObj <- .clear_after(scmpObj, "counts")
    validObject(scmpObj)

    # return
    return(scmpObj)
}
