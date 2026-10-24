#' @title Select the significant terms of each gene
#'
#' @description
#' `sc.t.fit()` selects by stepwise regression, as maSigPro does, the
#' significant terms of the full model of each gene with a non-flat profile
#' ([sc.p.vector()]), and keeps the coefficients, p-values and t statistics
#' of the selected models.
#'
#' @param selection_method The stepwise method: `"backward"`, `"forward"`,
#'   `"two.ways.backward"` or `"two.ways.forward"`.
#' @param p_value The level at which a term enters or stays in a model; by
#'   default that of [sc.p.vector()].
#' @param nvar_correction Whether to divide `p_value` by the number of terms
#'   of the design.
#' @param family The distribution of the error term, a [stats::family()]
#'   object; by default the family of [sc.p.vector()], with its link. When
#'   [sc.p.vector()] estimated one theta per gene (`theta_estimate = "gene"`),
#'   each gene is fitted with its own theta.
#' @param link The link of `family`, by default that of [sc.p.vector()]; when
#'   given, it must be the link of `family`. Choose the link in
#'   [sc.p.vector()].
#' @param offset Whether to add the offset of each bin, the logarithm of its
#'   size factor, which needs the log link. The fits use the offsets that
#'   [sc.p.vector()] stored in the `Design` slot; when it stored none
#'   (`offset = FALSE`), they are computed here in the same way.
#'
#' @importFrom BiocParallel bplapply bpnworkers SerialParam
#' @importFrom stats influence.measures
#'
#' @return The object with the selected model of each gene with a non-flat
#' profile in the `Estimate` slot (class [Estimates-class]), read by
#' [showSol()] and the other result readers.
#'
#' @seealso [showSol()] and the other readers of the results;
#' [Estimates-class].
#'
#' @references
#' Conesa A, Nueda MJ, Ferrer A, Talón M (2006). maSigPro: a method to
#' identify significantly differential expression profiles in time-course
#' microarray experiments. Bioinformatics 22(9), 1096-1102.
#' \doi{10.1093/bioinformatics/btl056}
#'
#' Nueda MJ, Tarazona S, Conesa A (2014). Next maSigPro: updating maSigPro
#' bioconductor package for RNA-seq time series. Bioinformatics 30(18),
#' 2598-2602. \doi{10.1093/bioinformatics/btu333}
#'
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#'
#' @keywords regression models
#' @family workflow steps
#' @inheritParams scmp-params
#' @inheritParams sc.p.vector
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' # Keep the terms with a p-value below 0.01
#' scmp <- sc.t.fit(scmp.ob, p_value = 0.01, verbose = FALSE)
#' head(showSol(scmp)[, c("p-value", "R-squared")])
#'
#' @export
sc.t.fit <- function(scmpObj,
                     selection_method = "backward",
                     p_value = scmpObj@Parameters@p_value,
                     nvar_correction = FALSE,
                     family = scmpObj@Parameters@distribution,
                     epsilon = scmpObj@Parameters@epsilon,
                     offset = scmpObj@Parameters@offset,
                     verbose = TRUE,
                     max_it = scmpObj@Parameters@max_it,
                     link = scmpObj@Parameters@link,
                     BPPARAM = SerialParam()) {
    # Check Object Validity
    .check_scmp(scmpObj)

    # The genes come from sc.p.vector() and the counts from pb_counts()
    .check(length(scmpObj@Profile@non_flat) > 0,
        msg = paste(
            "No gene with a non-flat profile found. Please run 'sc.p.vector()'",
            "first; if it found none, try changing its parameters."
        )
    )
    .check_dense_counts(scmpObj,
        bins = rownames(scmpObj@Design@predictor_matrix)
    )

    # Check the link when it is given (its default is the stored link)
    link_given <- !missing(link)
    if (link_given) {
        .check_choice(link, c("log", "identity"),
            msg = "link function should be either 'log' or 'identity'"
        )
    }
    selection_method <- .check_choice(selection_method,
        c("backward", "forward", "two.ways.backward", "two.ways.forward"),
        msg = paste(
            "stepwise method must be one of backward, forward,",
            "two.ways.backward, two.ways.forward"
        )
    )

    # The family keeps its own link (see .fit_family())
    family <- .fit_family(family, link, link_given, offset)
    link <- family$link

    # Transfer Data
    dis <- scmpObj@Design@predictor_matrix
    groups.vector <- scmpObj@Design@groups.vector
    groups.vector <- c(
        groups.vector[nchar(groups.vector) == min(nchar(groups.vector))][1],
        groups.vector
    )
    alloc <- scmpObj@Design@assignment_matrix
    G <- scmpObj@Parameters@g

    dat.all <- assay(scmpObj@Dense, "bulk.counts", withDimnames = FALSE)
    dat <- dat.all[scmpObj@Profile@non_flat, , drop = FALSE]
    dat <- rbind(c(rep(1, ncol(dat))), dat)
    dat <- dat[, as.character(rownames(dis))]

    g <- (dim(dat)[1] - 1)
    n <- dim(dat)[2]
    p <- dim(dis)[2]
    vars.in <- colnames(dis)
    sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
    influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
    rownames(influ.info) <- rownames(dis)

    # Why: a Bonferroni correction for the number of terms tested in
    # each gene.
    if (nvar_correction) {
        p_value <- p_value / ncol(dis)
    }

    # The offset of each bin: the one that sc.p.vector() computed from every
    # gene it fitted, or the same computed here when it stored none
    if (offset) {
        offsetData <- scmpObj@Design@offset
        if (length(offsetData) == 0L) {
            bins <- as.character(rownames(dis))
            counts <- as.matrix(dat.all[, bins, drop = FALSE])
            offsetData <- .bin_offset(
                .fit_genes(counts, scmpObj@Parameters@min_na)
            )
            scmpObj@Design@offset <- offsetData
        }
    } else {
        offsetData <- NULL
        scmpObj@Design@offset <- numeric(0)
    }

    # Collect all the y
    y_input <- lapply(2:(g + 1), function(i) as.numeric(dat[i, ]))
    names(y_input) <- rownames(dat)[-1]

    # Select the covariates
    if (selection_method == "backward") {
        step_fun <- sc.stepback
    } else if (selection_method == "forward") {
        step_fun <- sc.stepfor
    } else if (selection_method == "two.ways.backward") {
        step_fun <- sc.two.ways.stepback
    } else if (selection_method == "two.ways.forward") {
        step_fun <- sc.two.ways.stepfor
    }

    if (verbose) {
        message(
            "Selecting the terms of ", g, " genes with ",
            class(BPPARAM)[1L], " (workers: ", bpnworkers(BPPARAM), ")"
        )
    }

    # The fits name the terms by their position (see .fit_frame()); vars.in
    # keeps the design names for the columns of the results
    fit_dis <- .fit_frame(dis)
    theta_gene <- NULL
    if (identical(family$theta_estimate, "gene")) {
        theta_gene <- family$theta_gene[names(y_input)]
    }
    result_list <- bplapply(names(y_input), .tfit_gene,
        y_input = y_input, step_fun = step_fun, dis = fit_dis,
        family = family, epsilon = epsilon, offset_data = offsetData,
        vars_in = names(fit_dis), alfa = p_value, influ_info = influ.info,
        max_it = max_it, theta = theta_gene, BPPARAM = BPPARAM
    )

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

    # A gene whose selected model has no significant coefficient has no row
    # in sol, and none in the coefficient and t-score matrices (as maSigPro)
    in_sol <- !vapply(sol.list, is.null, logical(1))
    feature_names <- feature_names[in_sol]
    coeff.list <- coeff.list[in_sol]
    t.list <- t.list[in_sol]
    # Get influ.info
    influ.info.list <- lapply(result_list, function(element) {
        return(element[["influ.info"]])
    })

    # Keep the genes that have influential observations
    influ.info.list <- influ.info.list[!vapply(
        influ.info.list, is.logical,
        logical(1)
    )]

    # Lapply to remove column 1
    influ.info.list <- lapply(influ.info.list, function(element) {
        return(element[, -1, drop = FALSE])
    })

    # Create scmpObjframe
    sol <- do.call("rbind", sol.list)
    .check(!is.null(sol),
        msg = paste0(
            "No gene has a coefficient significant at p_value = ", p_value,
            ". Try a larger 'p_value'."
        )
    )
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
        if (!is.null(groups.vector) && !is.null(alloc)) {
            groups <- colnames(alloc)[3:ncol(alloc)]
            degree <- (length(groups.vector) / length(groups)) -
                1
            for (w in seq_len(nrow(coefficients))) {
                A <- NULL
                col.names <- NULL
                for (l in seq_along(groups)) {
                    B <- .path_coeffs(
                        coefficients = coefficients[w, ],
                        groups.vector = groups.vector, group = groups[l]
                    )
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
            message(
                "\nInfluence: ", ncol(influ.info), " genes with influential ",
                "observations detected. Model validation for these genes is ",
                "recommended"
            )
        }
    } else {
        if (verbose) {
            message("\nNo genes with influential observations detected")
        }
        influ.info <- matrix(data = NA, nrow = 0, ncol = 0)
    }

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
    scmpObj@Parameters@max_it <- as.integer(max_it)
    scmpObj@Parameters@offset <- offset
    scmpObj@Parameters@epsilon <- epsilon
    scmpObj@Parameters@selection_method <- selection_method
    scmpObj@Parameters@distribution <- family
    scmpObj@Parameters@link <- link
    validObject(scmpObj)

    return(scmpObj)
}

# Stepwise term selection for one gene with step_fun (one of the four maSigPro
# selection functions), the full and intercept-only fits, and the summaries
# of extract_fitting(). With `theta` (one per gene) the gene's fits use
# theta[[gene_name]]. sc.t.fit() runs it through bplapply(), so it takes all
# its data as arguments: SnowParam workers are new R sessions.
.tfit_gene <- function(gene_name, y_input, step_fun, dis, family, epsilon,
                       offset_data, vars_in, alfa, influ_info, max_it,
                       theta = NULL) {
    y <- y_input[[gene_name]]
    if (!is.null(theta)) {
        family <- .nb_with_theta(family, theta[[gene_name]])
    }

    reg <- step_fun(
        y = y, d = as.data.frame(dis),
        alfa = alfa, family = family,
        epsilon = epsilon,
        useOffset = offset_data,
        useWeight = NULL,
        max_it = max_it
    )
    lmf <- glm(y ~ .,
        data = as.data.frame(dis),
        family = family, epsilon = epsilon,
        offset = offset_data,
        maxit = max_it
    )
    model.glm.0 <- glm(y ~ 1,
        family = family,
        epsilon = epsilon,
        offset = offset_data,
        maxit = max_it
    )
    return(extract_fitting(
        reg = reg, lmf = lmf,
        model.glm.0 = model.glm.0, dis = dis,
        family = family, name = gene_name,
        vars.in = vars_in, alfa = alfa,
        influ.info = influ_info
    ))
}
