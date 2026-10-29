#' @title Test the full model of each gene
#'
#' @description
#' `sc.p.vector()` fits the full polynomial model of each gene to its
#' pseudo-bulk counts and tests it against the intercept-only model. The
#' genes whose adjusted p-value is at most `p_value` have a non-flat profile;
#' [sc.t.fit()] selects the terms of their models.
#'
#' @param p_value The level of the adjusted p-values: the genes at or below it
#'   have a non-flat profile. It is also the default level of [sc.t.fit()].
#' @param mt_correction The method of [stats::p.adjust()] that adjusts the
#'   p-values for multiple testing; `"BH"` is that of Benjamini and Hochberg
#'   (1995).
#' @param min_na The minimum number of bins with a value (not `NA`) that a
#'   gene needs to be fitted, as `min.obs` of `p.vector()` in maSigPro. Only
#'   counts with missing values, such as microarray data given to
#'   [create_scmp()] with `use_as_bin = TRUE`, have fewer; it can be at most
#'   the number of bins of the design.
#' @param family The distribution of the error term, a [stats::family()]
#'   object: by default the negative binomial with `theta = 10` and the link
#'   `link`. A family given here keeps its own link.
#' @param link The link of the default family: `"log"`, the model of the
#'   paper, or `"identity"`. When `family` is given, `link` must be its link.
#'   The identity link needs `offset = FALSE`; under it the fitted mean of a
#'   gene can become negative, and [stats::glm()] then warns.
#' @param epsilon The convergence tolerance of [stats::glm()].
#' @param offset Whether to add the offset of each bin, the logarithm of its
#'   size factor (see Details), to the linear predictor; it needs the log link.
#' @param max_it The maximum number of iterations of [stats::glm()].
#'
#' @details
#' With `offset = TRUE`, the fits include an offset as Srivastava et al.
#' (2024), Section 2.3 and Supplementary Section S1.3, describe it:
#' "scMaSigPro incorporates an offset value to the Poly-GLM model to adjust
#' for disparities in total read counts per bin. This offset is calculated as
#' the logarithm of size factors for each bin". The size factors are those of
#' the median-of-ratios method (Anders and Huber 2010; Love et al. 2014) on
#' the raw counts: "The size factor for each bin is the median of the ratios
#' of its counts to the geometric means for each gene". Only the genes with a
#' count above 0 in every bin enter the medians, and the fit stops when there
#' is none: fewer, larger bins ([sc.squeeze()]) or `offset = FALSE` avoid it.
#' The genes are fitted on their raw counts, without a pseudo-count:
#' "scMaSigPro can directly process raw counts". The offsets are stored in
#' the `Design` slot, where [sc.t.fit()] finds them.
#'
#' Genes are left out of the fit when they have fewer than `min_na` values,
#' when all their counts are zero, and when they are below `min_counts` or
#' `min_bins_expressed`. They have no p-value and do not count in the
#' multiple-testing correction; [scmpResults()] gives them `NA`. The size
#' factors are computed before `min_counts` and `min_bins_expressed` are
#' applied, so these two filters change the results of the genes they keep
#' only through the correction.
#'
#' The test of the full model is the first step of the method of Srivastava
#' et al. (2024), Section 2.4: "scMaSigPro implements a two-step strategy
#' similar to the maSigPro method", in which the full model of each gene is
#' "evaluated against a simpler intercept-only model via hypothesis testing
#' based on the log-likelihood ratio statistic". The default,
#' `test = "chisq"`, is that test as maSigPro computes it for families other
#' than the Gaussian: the statistic is divided by the estimated dispersion and
#' compared with the chi-squared distribution, as if the dispersion were
#' known. With few bins, and so few residual degrees of freedom, the estimate
#' varies from gene to gene and the test calls more genes significant than
#' its level. `test = "F"`, the quasi-likelihood F test, accounts for the
#' estimation of the dispersion through the residual degrees of freedom of
#' the F distribution.
#'
#' @param min_counts Minimum total count of a gene, summed over the bins of
#' the design. Genes with fewer counts are left out of the fit. The default,
#' 0, keeps every gene.
#' @param min_bins_expressed Minimum number of bins of the design in which a
#' gene has a count above 0. Genes expressed in fewer bins are left out of
#' the fit. The default, 0, keeps every gene; it can be at most the number of
#' bins.
#' @param theta_estimate For a negative binomial `family`: \code{"none"}
#' (the default) fits every gene with the theta of `family` (10 by default,
#' as maSigPro). \code{"gene"} estimates theta by maximum likelihood on each
#' gene's full model and fits the gene with it; \code{"global"} fits every
#' gene with the median of these estimates. The estimates are kept in the
#' `distribution` parameter, for [sc.t.fit()] and [plotDiagnostics()]:
#' `theta_gene` holds one per gene and `theta_estimate` the choice.
#' @param test The test of each gene's full model against its intercept-only
#' model: \code{"chisq"} (the default) or \code{"F"}. Both divide the
#' difference in deviance between the two models by the dispersion of the
#' full model, which, for the negative binomial and Gaussian families, is
#' estimated as its Pearson statistic divided by its residual degrees of
#' freedom. \code{"chisq"} compares the result with the chi-squared
#' distribution, as maSigPro does. \code{"F"} also divides it by the
#' difference in degrees of freedom and compares it with the F distribution
#' on the residual degrees of freedom of the full model: the F test of
#' [stats::anova.glm()]. With the Gaussian family both values give the F
#' test, the test of maSigPro for that family. \code{"F"} needs a family
#' whose dispersion is estimated, so not [stats::binomial()] or
#' [stats::poisson()]. The value is kept in the `test` parameter. See
#' Details.
#'
#' @importFrom stats anova dist glm make.link median p.adjust
#' @importFrom BiocParallel bplapply bpnworkers SerialParam
#' @importFrom MASS negative.binomial
#'
#' @return The object with the p-values and the genes with a non-flat
#' profile in the `Profile` slot (class [VariableProfiles-class]), the
#' offsets in the `Design` slot and the settings of the fit in the
#' `Parameters` slot. The `Estimate` and `Significant` slots,
#' computed from earlier tests, are emptied.
#'
#' @seealso [scmpResults()] for the p-values; [VariableProfiles-class].
#'
#' @references
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' Conesa A, Nueda MJ, Ferrer A, Talón M (2006). maSigPro: a method to
#' identify significantly differential expression profiles in time-course
#' microarray experiments. Bioinformatics 22(9), 1096-1102.
#' \doi{10.1093/bioinformatics/btl056}
#'
#' Nueda MJ, Tarazona S, Conesa A (2014). Next maSigPro: updating maSigPro
#' bioconductor package for RNA-seq time series. Bioinformatics 30(18),
#' 2598-2602. \doi{10.1093/bioinformatics/btu333}
#'
#' Anders S, Huber W (2010). Differential expression analysis for sequence
#' count data. Genome Biology 11(10), R106.
#' \doi{10.1186/gb-2010-11-10-r106}
#'
#' Love MI, Huber W, Anders S (2014). Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. Genome Biology 15(12), 550.
#' \doi{10.1186/s13059-014-0550-8}
#'
#' Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing. Journal of the Royal
#' Statistical Society: Series B 57(1), 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Venables WN, Ripley BD (2002). Modern Applied Statistics with S, fourth
#' edition. Springer, New York. ISBN 0-387-95457-0.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#'
#' @keywords regression models
#' @family workflow steps
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' # Test the full model of every gene with the F test
#' scmp <- sc.p.vector(scmp.ob, test = "F", verbose = FALSE)
#' sum(scmpResults(scmp)$non_flat)
#'
#' @export
sc.p.vector <- function(scmpObj, p_value = 0.05, mt_correction = "BH",
                        min_na = 6,
                        min_counts = 0,
                        min_bins_expressed = 0,
                        family = negative.binomial(theta = 10, link = link),
                        epsilon = 1e-8,
                        verbose = TRUE,
                        offset = TRUE,
                        max_it = 100,
                        link = "log",
                        theta_estimate = "none",
                        test = "chisq",
                        BPPARAM = SerialParam()) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # The design comes from sc.set.poly() and the counts from pb_counts()
    .check(ncol(scmpObj@Design@predictor_matrix) > 0,
        msg = "No polynomial design found. Please run 'sc.set.poly()', first."
    )
    .check_dense_counts(scmpObj,
        bins = rownames(scmpObj@Design@predictor_matrix)
    )

    # Extract from s4
    dis <- .fit_frame(scmpObj@Design@predictor_matrix)
    groups.vector <- scmpObj@Design@groups.vector
    alloc <- scmpObj@Design@assignment_matrix

    # Convert 'scmpObj' to matrix and select relevant columns based on 'design' rows
    dat <- as.matrix(assay(scmpObj@Dense, "bulk.counts", withDimnames = FALSE))
    dat <- dat[, as.character(rownames(dis)), drop = FALSE]
    G <- nrow(dat)

    # The link of the default family, checked before the family is built
    link_given <- !missing(link)
    .check_choice(link, c("log", "identity"),
        msg = "link function should be either 'log' or 'identity'"
    )

    # The fits use the link functions of the family
    family <- .fit_family(family, link, link_given, offset)
    link <- family$link
    theta_estimate <- .check_choice(theta_estimate,
        c("none", "gene", "global"),
        msg = "'theta_estimate' must be one of 'none', 'gene' or 'global'"
    )
    .check(theta_estimate == "none" || .is_nb_family(family),
        msg = paste(
            "'theta_estimate' needs a negative binomial family, such as",
            "MASS::negative.binomial(10)."
        )
    )
    test <- .check_choice(test, c("chisq", "F"),
        msg = "'test' must be either 'chisq' or 'F'"
    )
    .check(test == "chisq" || !.fixed_dispersion(family),
        msg = paste0(
            "test = 'F' needs a family whose dispersion is estimated: the '",
            family$family, "' family fixes it at 1. Use test = 'chisq'."
        )
    )

    # Add check
    .check(min_na <= ncol(dat),
        msg = paste(
            min_na,
            "for 'min_na' is too high. Try lowering the threshold."
        )
    )
    .check(
        is.numeric(min_counts) && length(min_counts) == 1L &&
            !is.na(min_counts) && min_counts >= 0,
        msg = "'min_counts' must be a single number of at least 0."
    )
    .check(
        is.numeric(min_bins_expressed) && length(min_bins_expressed) == 1L &&
            !is.na(min_bins_expressed) && min_bins_expressed >= 0 &&
            min_bins_expressed <= ncol(dat),
        msg = paste0(
            "'min_bins_expressed' must be a single number from 0 to ",
            ncol(dat), ", the number of bins of the design."
        )
    )

    # Removing rows with many missings and rows with all zeros
    dat <- .fit_genes(dat, min_na)

    # The count filter. It is applied after the size factors
    # below, so it changes the results of the genes it keeps only through
    # the multiple-testing correction; a threshold of 0 keeps every gene
    keep <- rep(TRUE, nrow(dat))
    if (min_counts > 0) {
        keep <- keep & rowSums(dat, na.rm = TRUE) >= min_counts
    }
    if (min_bins_expressed > 0) {
        keep <- keep & rowSums(dat > 0, na.rm = TRUE) >= min_bins_expressed
    }

    # The offset of each bin: the logarithm of its size factor (Srivastava
    # et al. (2024), Section 2.3), computed once here for sc.t.fit() too
    # Why: under the log link, log(s) in the linear predictor multiplies
    # the mean of a bin by its size factor s, so the genes are compared at
    # equal depth; sc.t.fit() reuses the same offsets.
    if (offset) {
        offsetData <- .bin_offset(dat)
        scmpObj@Design@offset <- offsetData
    } else {
        offsetData <- NULL
        scmpObj@Design@offset <- numeric(0)
    }

    # The genes below min_counts or min_bins_expressed are not fitted
    if (verbose && !all(keep)) {
        message(
            "Leaving out ", sum(!keep), " of ", length(keep), " genes ",
            "below 'min_counts' or 'min_bins_expressed'"
        )
    }
    dat <- dat[keep, , drop = FALSE]

    # Get dimensions for the input
    g <- dim(dat)[1]
    n <- dim(dat)[2]
    p <- dim(dis)[2]
    sc.p.vector <- vector(mode = "numeric", length = g)

    if (verbose) {
        message(
            "Fitting the full model to ", g, " genes with ",
            class(BPPARAM)[1L], " (workers: ", bpnworkers(BPPARAM), ")"
        )
    }

    # Estimated theta: one per gene on its full model, then
    # each gene's own (gene) or their median (global)
    theta_gene <- NULL
    if (theta_estimate != "none") {
        est <- bplapply(seq_len(g), .theta_gene,
            dat = dat, dis = dis, family = family, epsilon = epsilon,
            offset_data = offsetData, max_it = max_it, BPPARAM = BPPARAM
        )
        est <- matrix(as.numeric(unlist(est)),
            ncol = 4L, byrow = TRUE,
            dimnames = list(
                rownames(dat), c("theta", "at_bound", "converged", "n_alt")
            )
        )
        thetas <- est[, "theta"]
        names(thetas) <- rownames(dat)
        if (verbose && g > 0L) {
            message(
                "Estimated theta (", theta_estimate, "): median ",
                signif(median(thetas), 4), ", range ",
                signif(min(thetas), 4), " to ", signif(max(thetas), 4),
                "; ", sum(est[, "at_bound"] == 1), " of ", g,
                " at a bound, ", sum(est[, "converged"] == 0),
                " not converged"
            )
        }
        if (theta_estimate == "gene") {
            theta_gene <- thetas
        } else if (g > 0L) {
            family <- .nb_with_theta(family, median(thetas))
        }
        family$theta_estimate <- theta_estimate
        family$theta_gene <- thetas
    }

    p.vector.list <- bplapply(seq_len(g), .pvector_gene,
        dat = dat, dis = dis, family = family, epsilon = epsilon,
        offset_data = offsetData, max_it = max_it, theta = theta_gene,
        f_test = test == "F", BPPARAM = BPPARAM
    )

    names(p.vector.list) <- rownames(dat)
    sc.p.vector <- unlist(p.vector.list, recursive = TRUE, use.names = TRUE)
    if (is.null(sc.p.vector)) {
        # No gene is left after the filters above
        sc.p.vector <- numeric(0)
    }
    #----------------------------------------------------------------------
    # Correct p-values using FDR correction and select significant genes
    # Why: every fitted gene is tested; the correction (Benjamini and
    # Hochberg by default) runs over the fitted genes only.
    p.adjusted <- unlist(
        p.adjust(sc.p.vector,
            method = mt_correction,
            n = length(sc.p.vector)
        ),
        recursive = TRUE, use.names = TRUE
    )
    names(p.adjusted) <- names(sc.p.vector)
    genes.selected <- rownames(dat)[which(p.adjusted <= p_value)]
    # Why: the correction keeps the order of the p-values, so the selected
    # genes have the smallest ones and this is the largest of them.
    FDR <- sort(sc.p.vector)[length(genes.selected)]

    # Subset the expression values of significant genes
    SELEC <- dat[rownames(dat) %in% genes.selected, , drop = FALSE]

    # Prepare 'sc.p.vector' for output
    names(sc.p.vector) <- rownames(dat)

    # Add Data to the class. With no significant gene the profile has no
    # non-flat gene, and the profile and parameters are written all the same
    profile.obj <- new("VariableProfiles",
        non_flat = as.character(rownames(SELEC)),
        p_values = sc.p.vector,
        adj_p_values = p.adjusted,
        fdr = FDR
    )

    # Update Slot
    scmpObj@Profile <- profile.obj
    scmpObj <- .clear_after(scmpObj, "tests")

    # Update Parameter Slot
    scmpObj@Parameters@max_it <- as.integer(max_it)
    scmpObj@Parameters@offset <- offset
    scmpObj@Parameters@p_value <- p_value
    scmpObj@Parameters@min_na <- min_na
    scmpObj@Parameters@min_counts <- as.numeric(min_counts)
    scmpObj@Parameters@min_bins_expressed <- as.numeric(min_bins_expressed)
    scmpObj@Parameters@g <- g
    scmpObj@Parameters@mt_correction <- mt_correction
    scmpObj@Parameters@test <- test
    scmpObj@Parameters@epsilon <- epsilon
    scmpObj@Parameters@distribution <- family
    scmpObj@Parameters@link <- link

    if (nrow(SELEC) == 0) {
        message("No significant genes detected. Try changing parameters.")
    }
    validObject(scmpObj)

    return(scmpObj)
}

# The family of the fits of sc.p.vector() and sc.t.fit(). glm() uses the
# link functions of a family (linkfun, linkinv, mu.eta), so a family is used
# with its own link, and renaming its `link` entry would only relabel it.
# A relabelled family, such as the one an object fitted with
# link = "identity" by scMaSigPro 0.0.5 or earlier holds, is refused.
# `link`, when given, must be the link of the family. The size-factor
# offset is added to the linear predictor, where it scales the mean under
# the log link only, so offset = TRUE needs the log link.
.fit_family <- function(family, link, link_given, offset) {
    .check(inherits(family, "family"),
        msg = paste(
            "'family' must be a family object, such as",
            "MASS::negative.binomial(theta = 10)."
        )
    )
    .check(!.relabelled(family), msg = .relabelled_message(family))
    .check(!link_given || identical(family$link, link),
        msg = paste0(
            "'family' has the '", family$link, "' link, not link = '", link,
            "'. To fit with the '", link, "' link, give a family built ",
            "with it, such as MASS::negative.binomial(theta = 10, link = '",
            link, "'), or leave 'family' out of sc.p.vector()."
        )
    )
    .check(!isTRUE(offset) || identical(family$link, "log"),
        msg = paste0(
            "The size-factor offset applies to the log link only: with the '",
            family$link, "' link, pass offset = FALSE."
        )
    )
    family
}

# TRUE for a family whose dispersion summary.glm() and anova() fix instead
# of estimating: one that declares it, or the binomial and Poisson families
.fixed_dispersion <- function(family) {
    dispersion <- family[["dispersion"]]
    (!is.null(dispersion) && !is.na(dispersion)) ||
        family$family %in% c("binomial", "poisson")
}

# The link whose functions `family` has: its own link label when the
# functions are those of make.link() for that label, otherwise the first of
# the links of make.link() whose linkfun and linkinv give the same values on
# a few test points, or NA when there is none (a link make.link() does not
# know, such as a power link, is not checked).
.true_link <- function(family) {
    x <- c(0.25, 0.5, 0.75)
    links <- c(
        family$link, "log", "identity", "logit", "probit", "cauchit",
        "cloglog", "sqrt", "inverse", "1/mu^2"
    )
    for (name in unique(links)) {
        same <- tryCatch(
            {
                ref <- make.link(name)
                isTRUE(all.equal(family$linkfun(x), ref$linkfun(x))) &&
                    isTRUE(all.equal(family$linkinv(x), ref$linkinv(x)))
            },
            error = function(e) FALSE,
            warning = function(w) FALSE
        )
        if (same) {
            return(name)
        }
    }
    NA_character_
}

# TRUE for a family whose link functions are those of another link than its
# label: earlier versions renamed the link of a family without changing its
# functions, so an object fitted by them with a link other than
# the family's holds such a family
.relabelled <- function(family) {
    if (!inherits(family, "family")) {
        return(FALSE)
    }
    actual <- .true_link(family)
    !is.na(actual) && !identical(actual, family$link)
}

.relabelled_message <- function(family) {
    paste0(
        "'family' is labelled with the '", family$link, "' link but has the ",
        "link functions of the '", .true_link(family), "' link. scMaSigPro ",
        "0.0.5 and earlier renamed the link of the family without changing ",
        "its functions, so an object fitted by them with link = \"",
        family$link, "\" holds fits made with the '", .true_link(family),
        "' link: run sc.p.vector() on it again."
    )
}

# The genes that sc.p.vector() fits, before its count filter: the rows of
# dat with at least min_na values (only counts with NA values have fewer)
# and a count above 0. sc.t.fit() computes the offset from the same genes
# when sc.p.vector() stored none.
.fit_genes <- function(dat, min_na) {
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min_na, , drop = FALSE]
    sumatot <- apply(dat, 1, sum)
    counts0 <- which(sumatot == 0)
    if (length(counts0) > 0) {
        dat <- dat[-counts0, , drop = FALSE]
    }
    dat
}

# The design as a data frame whose columns are named by position: x1, x2, ...
# The per-gene fits of sc.p.vector() and sc.t.fit() use these names in their
# formulas (y ~ .), so any design names work, including names with spaces or
# `|`, repeated names and "y"; extract_fitting() then matches the terms of a
# fit to the columns of the design by position.
.fit_frame <- function(dis) {
    dis <- as.data.frame(dis)
    names(dis) <- paste0("x", seq_len(ncol(dis)))
    dis
}

# Full-model p-value of gene i (row i of dat): glm(y ~ .) against glm(y ~ 1),
# with an F test for the Gaussian family and, for the other families, a
# chi-squared test, or an F test when f_test is TRUE. With `theta` (one
# estimate per gene) both models of gene i use theta[[i]]. sc.p.vector()
# runs it through bplapply(), so it takes all its data as arguments:
# SnowParam workers are new R sessions.
.pvector_gene <- function(i, dat, dis, family, epsilon, offset_data, max_it,
                          theta = NULL, f_test = FALSE) {
    y <- as.numeric(dat[i, ])
    if (!is.null(theta)) {
        family <- .nb_with_theta(family, theta[[i]])
    }

    # Set full Model
    model.glm <- glm(y ~ .,
        data = dis, family = family, epsilon = epsilon,
        offset = offset_data,
        maxit = max_it
    )
    sc_p_val <- NA
    # Why: when the intercept-only model fits exactly, the full model
    # explains nothing more.
    if (model.glm$null.deviance == 0) {
        sc_p_val <- 1
    } else {
        model.glm.0 <- glm(y ~ 1,
            family = family, epsilon = epsilon,
            offset = offset_data,
            maxit = max_it
        )

        # Perform ANOVA or Chi-square test based on the distribution
        if (family$family == "gaussian") {
            # Why: the Gaussian family estimates its variance, and the F test
            # takes that estimate into account.
            test <- anova(model.glm.0, model.glm, test = "F")
            sc_p_val <- ifelse(is.na(test[6][2, 1]), 1, test[6][2, 1])
        } else if (f_test) {
            # Without residual degrees of freedom anova() gives no F
            # statistic: the p-value is then 1, as with the chi-squared test
            test <- anova(model.glm.0, model.glm, test = "F")
            p <- if ("Pr(>F)" %in% names(test)) test[2L, "Pr(>F)"] else NA
            sc_p_val <- ifelse(is.na(p), 1, p)
        } else {
            # Why: the log-likelihood ratio test of Section 2.4 of the paper,
            # with the deviance difference scaled by the dispersion of the full
            # model.
            test <- anova(model.glm.0, model.glm, test = "Chisq")
            sc_p_val <- ifelse(is.na(test[5][2, 1]), 1, test[5][2, 1])
        }
    }

    return(sc_p_val)
}
