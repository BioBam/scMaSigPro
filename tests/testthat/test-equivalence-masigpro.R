# Oracle tests: with no offset, scMaSigPro must reproduce maSigPro exactly.
# test-reproduce-masigpro.R covers the Gaussian/identity case on the full
# data.abiotic with backward selection. This file extends the oracle to every
# selection method and polynomial degree (on a 200-gene subset, for speed), and
# to the negative binomial family on maSigPro's NBdata.

quiet <- function(expr) {
    invisible(utils::capture.output(res <- expr))
    res
}

group_column <- function(design, group_cols) {
    md <- as.data.frame(design)
    md$Group <- apply(md[, group_cols, drop = FALSE], 1, function(x) names(x)[x == 1])
    md[, !(colnames(md) %in% group_cols), drop = FALSE]
}

methods <- c("backward", "forward", "two.ways.backward", "two.ways.forward")

# Every method runs with the default BPPARAM (SerialParam()). Before
# BiocParallel, the progress-bar code of the three non-backward methods
# errored on R >= 4.3 (ledger T1).
fit_sc <- function(ob, method) {
    sc.t.fit(ob, selection_method = method, verbose = FALSE)
}

expect_tfit_matches <- function(sc, tstep, non_flat) {
    keep <- rownames(tstep$sol)
    expect_identical(showSol(sc), tstep$sol)
    expect_identical(showCoeff(sc)[keep, , drop = FALSE], tstep$coefficients)
    expect_identical(as.data.frame(showGroupCoeff(sc))[keep, , drop = FALSE], tstep$group.coeffs)
    expect_identical(as.data.frame(showTS(sc))[keep, , drop = FALSE], tstep$t.score)
    # Known divergence (ledger T2): scMaSigPro keeps coefficient and t-score
    # rows for genes whose selected model has no significant term; maSigPro
    # drops them from every output.
    expect_identical(setdiff(rownames(showCoeff(sc)), keep), setdiff(non_flat, keep))
    infl_sc <- sc@Estimate@influential
    infl_ma <- tstep$influ.info
    # maSigPro drops influ.info to a named vector when exactly one gene is
    # influential (its `[, -1]` has no drop = FALSE); scMaSigPro keeps a matrix.
    if (is.null(dim(infl_ma))) {
        infl_ma <- matrix(infl_ma, ncol = 1, dimnames = list(names(infl_ma), colnames(infl_sc)))
    }
    expect_identical(ncol(infl_sc), ncol(infl_ma))
    if (ncol(infl_ma) > 0) {
        expect_identical(infl_sc, infl_ma)
    }
}

# ---- Gaussian, data.abiotic subset -------------------------------------------

abiotic_case <- function(degree) {
    fixture(paste0("oracle_abiotic_", degree), function() {
        e <- new.env()
        utils::data("data.abiotic", "edesign.abiotic", package = "maSigPro", envir = e)
        counts <- as.matrix(e$data.abiotic[seq_len(200), ])
        meta <- group_column(e$edesign.abiotic, c("Control", "Cold", "Heat", "Salt"))
        design <- maSigPro::make.design.matrix(e$edesign.abiotic, degree = degree)
        fit <- quiet(maSigPro::p.vector(counts, design, Q = 0.05, MT.adjust = "BH", min.obs = 20))
        ob <- create_scmp(
            counts = counts, cell_data = meta,
            ptime_col = "Time", path_col = "Group", use_as_bin = TRUE
        )
        ob <- sc.set.poly(ob, poly_degree = degree)
        ob <- sc.p.vector(ob,
            min_na = 20, verbose = FALSE, link = "identity", offset = FALSE,
            max_it = 25, epsilon = 0.00001, family = gaussian()
        )
        list(design = design, fit = fit, ob = ob)
    })
}

for (degree in 1:3) {
    test_that(paste0("Gaussian degree ", degree, ": design and p.vector match maSigPro"), {
        case <- abiotic_case(degree)
        expect_equal(predictors(case$ob), as.matrix(case$design$dis))
        expect_equal(pathAssign(case$ob), as.matrix(case$design$edesign))
        expect_identical(case$ob@Design@groups.vector, case$design$groups.vector)
        expect_identical(as.vector(case$ob@Profile@p_values), as.vector(case$fit$p.vector))
        expect_identical(as.vector(case$ob@Profile@adj_p_values), as.vector(case$fit$p.adjusted))
        expect_identical(case$ob@Profile@non_flat, rownames(case$fit$SELEC))
    })

    for (method in methods) {
        test_that(paste0("Gaussian degree ", degree, ", ", method, ": T.fit matches maSigPro"), {
            case <- abiotic_case(degree)
            tstep <- quiet(maSigPro::T.fit(case$fit, step.method = method, alfa = 0.05))
            expect_tfit_matches(fit_sc(case$ob, method), tstep, case$ob@Profile@non_flat)
        })
    }
}

# ---- Negative binomial, NBdata (theta = 10, log link, no offset) ------------

nb_case <- function() {
    fixture("oracle_nb", function() {
        e <- new.env()
        utils::data("NBdata", "NBdesign", package = "maSigPro", envir = e)
        design <- maSigPro::make.design.matrix(e$NBdesign, degree = 2)
        fit <- quiet(maSigPro::p.vector(e$NBdata, design,
            Q = 0.05, MT.adjust = "BH", min.obs = 6, counts = TRUE, theta = 10
        ))
        ob <- create_scmp(
            counts = as.matrix(e$NBdata), cell_data = group_column(e$NBdesign, c("Group.1", "Group.2")),
            ptime_col = "Time", path_col = "Group", use_as_bin = TRUE
        )
        ob <- sc.set.poly(ob, poly_degree = 2)
        ob <- sc.p.vector(ob,
            min_na = 6, verbose = FALSE, link = "log", offset = FALSE,
            family = MASS::negative.binomial(theta = 10), epsilon = 1e-5, max_it = 25
        )
        list(fit = fit, ob = ob)
    })
}

test_that("negative binomial p.vector matches maSigPro (counts = TRUE, theta = 10)", {
    case <- nb_case()
    expect_identical(as.vector(case$ob@Profile@p_values), as.vector(case$fit$p.vector))
    expect_identical(as.vector(case$ob@Profile@adj_p_values), as.vector(case$fit$p.adjusted))
    expect_identical(case$ob@Profile@non_flat, rownames(case$fit$SELEC))
    expect_length(case$ob@Profile@non_flat, 21L)
})

for (method in methods) {
    test_that(paste0("negative binomial, ", method, ": T.fit matches maSigPro"), {
        case <- nb_case()
        tstep <- quiet(maSigPro::T.fit(case$fit,
            step.method = method, alfa = 0.05,
            family = MASS::negative.binomial(theta = 10), epsilon = 1e-5
        ))
        expect_tfit_matches(fit_sc(case$ob, method), tstep, case$ob@Profile@non_flat)
    })
}

test_that("every selection method also runs with verbose = TRUE (ledger T1, q9)", {
    ob <- nb_case()$ob
    for (method in methods) {
        msgs <- capture_messages(fit <- sc.t.fit(ob, selection_method = method))
        expect_match(msgs[1], "with SerialParam (workers: 1)", fixed = TRUE)
        expect_identical(fit@Estimate, fit_sc(ob, method)@Estimate)
    }
})
