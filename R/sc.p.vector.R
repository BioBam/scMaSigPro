#' @title Perform fitting with full model.
#'
#' @description
#' Performs a regression fit for each gene taking all variables present in the
#' model.
#'
#' @importFrom stats anova dist glm median p.adjust
#' @importFrom BiocParallel bplapply bpnworkers SerialParam
#' @importFrom MASS negative.binomial
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param p_value Significance level used for variable selection in the stepwise
#' regression.
#' @param mt_correction A character string specifying the p-value correction
#' method.
#' @param min_na Minimum values needed per gene across cells to estimate the
#' model.
#' @param family Distribution of the error term, a \code{\link[stats]{family}}
#' object. The default is the negative binomial family with \code{theta = 10}
#' and the link given by \code{link}. A family given here is used with its
#' own link.
#' @param link Link function of the default family: "log" (the default, the
#' model of the paper) or "identity". When \code{family} is also given,
#' \code{link} must be the link of that family. The identity link needs
#' \code{offset = FALSE}; under it the fitted mean of a gene can become
#' negative, and \code{glm()} then warns.
#' @param epsilon Model convergence tolerance.
#' @param offset logical value specifying whether to use offset during fitting.
#' The offset, the size factor of each bin (or its logarithm, see
#' \code{log_offset}), is added to the linear predictor, where it scales the
#' fitted mean under the log link only: with any other link the fits stop
#' with an error unless \code{offset = FALSE}.
#' @param log_offset A logical value specifying whether to take the logarithm of
#' the offsets.
#' @param max_it Maximum number of iterations to fit the model.
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' @param BPPARAM A \code{BiocParallelParam} object that runs the per-gene
#' fits, such as
#' \code{\link[BiocParallel:MulticoreParam-class]{MulticoreParam}(2)} (Linux
#' and macOS) or \code{\link[BiocParallel:SnowParam-class]{SnowParam}(2)}
#' (any system, including Windows). The default,
#' \code{\link[BiocParallel:SerialParam-class]{SerialParam}()}, fits the genes
#' one after another. Create it with \code{progressbar = TRUE} for a progress
#' bar. The results do not depend on the backend.
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
                        family = negative.binomial(theta = 10, link = link),
                        epsilon = 1e-8,
                        verbose = TRUE,
                        offset = TRUE,
                        log_offset = FALSE,
                        max_it = 100,
                        link = "log",
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
    dat <- dat[, as.character(rownames(dis))]
    G <- nrow(dat)

    # The link of the default family, checked before the family is built
    link_given <- !missing(link)
    .check_choice(link, c("log", "identity"),
        msg = "link function should be either 'log' or 'identity'"
    )

    # The fits use the link functions of the family (ledger S6)
    family <- .fit_family(family, link, link_given, offset)
    link <- family$link

    # Add check
    .check(min_na <= ncol(dat),
        msg = paste(
            min_na,
            "for 'min_na' is too high. Try lowering the threshold."
        )
    )

    # Removing rows with many missings:
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min_na, , drop = FALSE]

    # Removing rows with all zeros:
    sumatot <- apply(dat, 1, sum)
    counts0 <- which(sumatot == 0)
    if (length(counts0) > 0) {
        dat <- dat[-counts0, , drop = FALSE]
    }

    # Get dimensions for the input
    g <- dim(dat)[1]
    n <- dim(dat)[2]
    p <- dim(dis)[2]
    sc.p.vector <- vector(mode = "numeric", length = g)

    # Calculate  offset
    if (offset) {
        dat <- dat + 1
        offsetData <- .size_factors(dat)
        if (log_offset) {
            offsetData <- log(offsetData)
        }
        scmpObj@Design@offset <- offsetData
    } else {
        offsetData <- NULL
    }

    if (verbose) {
        message(
            "Fitting the full model to ", g, " genes with ",
            class(BPPARAM)[1L], " (workers: ", bpnworkers(BPPARAM), ")"
        )
    }

    p.vector.list <- bplapply(seq_len(g), .pvector_gene,
        dat = dat, dis = dis, family = family, epsilon = epsilon,
        offset_data = offsetData, max_it = max_it, BPPARAM = BPPARAM
    )

    names(p.vector.list) <- rownames(dat)
    sc.p.vector <- unlist(p.vector.list, recursive = TRUE, use.names = TRUE)
    if (is.null(sc.p.vector)) {
        # No gene is left after the filters above
        sc.p.vector <- numeric(0)
    }
    #----------------------------------------------------------------------
    # Correct p-values using FDR correction and select significant genes
    p.adjusted <- unlist(
        p.adjust(sc.p.vector,
            method = mt_correction,
            n = length(sc.p.vector)
        ),
        recursive = TRUE, use.names = TRUE
    )
    names(p.adjusted) <- names(sc.p.vector)
    genes.selected <- rownames(dat)[which(p.adjusted <= p_value)]
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

    # Update Parameter Slot
    scmpObj@Parameters@log_offset <- log_offset
    scmpObj@Parameters@max_it <- as.integer(max_it)
    scmpObj@Parameters@offset <- offset
    scmpObj@Parameters@p_value <- p_value
    scmpObj@Parameters@min_na <- min_na
    scmpObj@Parameters@g <- g
    scmpObj@Parameters@mt_correction <- mt_correction
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
# with its own link, and renaming its `link` entry would only relabel it
# (ledger S6). `link`, when given, must be that link. The size-factor offset
# is added to the linear predictor, where it scales the mean under the log
# link only, so offset = TRUE needs the log link.
.fit_family <- function(family, link, link_given, offset) {
    .check(inherits(family, "family"),
        msg = paste(
            "'family' must be a family object, such as",
            "MASS::negative.binomial(theta = 10)."
        )
    )
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
# with an F test for the Gaussian family and a chi-squared test otherwise.
# sc.p.vector() runs it through bplapply(), so it takes all its data as
# arguments: SnowParam workers are new R sessions.
.pvector_gene <- function(i, dat, dis, family, epsilon, offset_data, max_it) {
    y <- as.numeric(dat[i, ])

    # Set full Model
    model.glm <- glm(y ~ .,
        data = dis, family = family, epsilon = epsilon,
        offset = offset_data,
        maxit = max_it
    )
    sc_p_val <- NA
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
            test <- anova(model.glm.0, model.glm, test = "F")
            sc_p_val <- ifelse(is.na(test[6][2, 1]), 1, test[6][2, 1])
        } else {
            test <- anova(model.glm.0, model.glm, test = "Chisq")
            sc_p_val <- ifelse(is.na(test[5][2, 1]), 1, test[5][2, 1])
        }
    }

    return(sc_p_val)
}
