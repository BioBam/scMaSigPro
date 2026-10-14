# The polynomial design and the path coefficients of R/poly_design.R on small
# designs worked out by hand, then against maSigPro's make.design.matrix() and
# reg.coeffs() on a grid of designs (skipped when maSigPro is not installed).

# A design as sc.set.poly() builds it: time, replicate index and one 0/1
# column per path, the paths in order of first appearance. `times` holds the
# time of every row of each path.
toy_design <- function(times, reps = 1L, interleave = FALSE) {
    path <- rep(rep(names(times), lengths(times)), each = reps)
    time <- rep(unlist(times, use.names = FALSE), each = reps)
    if (interleave) {
        ord <- order(ave(seq_along(path), path, FUN = seq_along))
        path <- path[ord]
        time <- time[ord]
    }
    cond <- paste(path, time)
    n <- length(path)
    m <- cbind(
        Time = time,
        Replicate = match(cond, unique(cond)),
        vapply(unique(path), function(p) as.numeric(path == p), numeric(n))
    )
    rownames(m) <- paste0(path, "_bin_", seq_along(path))
    m
}

two_paths <- toy_design(list(A = c(1, 2, 3), B = c(1, 2, 4)))

test_that(".poly_design() dummy-codes the paths and raises time", {
    d <- .poly_design(two_paths, degree = 2)
    expect_named(d, c("dis", "groups.vector", "edesign"))
    expect_identical(d$edesign, two_paths)
    expect_identical(
        colnames(d$dis), c("BvsA", "Time", "TimexB", "Time2", "Time2xB")
    )
    expect_identical(rownames(d$dis), rownames(two_paths))
    expect_identical(unname(d$dis[, "BvsA"]), c(0, 0, 0, 1, 1, 1))
    expect_identical(unname(d$dis[, "Time"]), c(1, 2, 3, 1, 2, 4))
    expect_identical(unname(d$dis[, "TimexB"]), c(0, 0, 0, 1, 2, 4))
    expect_identical(unname(d$dis[, "Time2"]), c(1, 4, 9, 1, 4, 16))
    expect_identical(unname(d$dis[, "Time2xB"]), c(0, 0, 0, 1, 4, 16))
    expect_identical(
        d$groups.vector, c("BvsA", "A", "BvsA", "A", "BvsA")
    )
})

test_that(".poly_design() names the time terms after the time column", {
    e <- two_paths
    colnames(e)[1] <- "scmp_binned_pseudotime"
    expect_identical(colnames(.poly_design(e, degree = 2)$dis), c(
        "BvsA", "scmp_binned_pseudotime", "scmp_binned_pseudotimexB",
        "scmp_binned_pseudotime2", "scmp_binned_pseudotime2xB"
    ))
})

test_that(".poly_design() with one path holds the powers of time only", {
    d <- .poly_design(toy_design(list(A = c(0, 1, 2, 3))), degree = 3)
    expect_identical(colnames(d$dis), c("Time", "Time2", "Time3"))
    expect_identical(unname(d$dis[, "Time3"]), c(0, 1, 8, 27))
    expect_identical(d$groups.vector, c("A", "A", "A"))
})

test_that(".poly_design() keeps an integer time for the first power", {
    e <- two_paths
    storage.mode(e) <- "integer"
    expect_type(.poly_design(e, degree = 1)$dis, "integer")
    d <- .poly_design(e, degree = 2)
    expect_type(d$dis, "double")
    expect_identical(unname(d$dis[, "Time"]), as.numeric(e[, "Time"]))
    e2 <- e
    rownames(e2) <- NULL
    expect_null(rownames(.poly_design(e2, degree = 2)$dis))
})

test_that(".poly_design() accepts a single earliest row", {
    # maSigPro's make.design.matrix() stops with "incorrect number of
    # dimensions" when only one row has the earliest time
    e <- toy_design(list(A = c(0, 1, 2), B = c(1, 2)))
    expect_identical(
        colnames(.poly_design(e, degree = 1)$dis),
        c("BvsA", "Time", "TimexB")
    )
})

test_that(".poly_design() needs a degree of at least 1", {
    for (degree in list(0, -1, 0.5, NA_real_, Inf, c(1, 2), "2")) {
        expect_error(.poly_design(two_paths, degree = degree),
            "'poly_degree' must be a number of at least 1",
            fixed = TRUE
        )
    }
    expect_error(sc.set.poly(splat_binned(), poly_degree = 0),
        "'poly_degree' must be a number of at least 1",
        fixed = TRUE
    )
    # A fractional degree uses its whole part, with any number of paths
    expect_identical(
        .poly_design(two_paths, degree = 2.5),
        .poly_design(two_paths, degree = 2L)
    )
    one_path <- toy_design(list(A = c(0, 1, 2, 3)))
    expect_identical(
        .poly_design(one_path, degree = 1.5),
        .poly_design(one_path, degree = 1)
    )
})

test_that(".path_coeffs() adds a path's terms to the reference ones", {
    gv <- .poly_design(two_paths, degree = 2)$groups.vector
    gv_fit <- c("A", gv)
    # beta0, BvsA, Time, TimexB, Time2, Time2xB
    cf <- c(10, 1, 2, 3, 4, 5)
    expect_identical(.path_coeffs(cf, gv_fit, "A"), c(10, 2, 4))
    expect_identical(.path_coeffs(cf, gv_fit, "B"), c(11, 5, 9))
    # Without the intercept's label (as plotTrend() passes it)
    expect_identical(.path_coeffs(cf, gv, "B"), c(11, 5, 9))
    # An unknown path gets the reference polynomial
    expect_identical(.path_coeffs(cf, gv_fit, "C"), c(10, 2, 4))
    # One row of the coefficient table, as sc.t.fit() passes it
    row <- as.data.frame(t(cf))
    colnames(row) <- c("beta0", paste0("beta", c(
        "BvsA", "Time", "TimexB", "Time2", "Time2xB"
    )))
    rownames(row) <- "gene1"
    expect_identical(.path_coeffs(row, gv_fit, "B"), c(11, 5, 9))
    cf[2] <- NA
    expect_identical(.path_coeffs(cf, gv_fit, "B"), c(NA, 5, 9))
    expect_identical(.path_coeffs(cf, gv_fit, "A"), c(10, 2, 4))
})

# ---- Oracle: maSigPro ------------------------------------------------------

oracle_designs <- function() {
    tricky <- toy_design(list(
        Y_15 = c(0.5, 1.25, 2.75), `Y_18_links_Y_15` = c(0.5, 3.1),
        `a b` = c(1, 2.2, 3.3, 4.4)
    ), reps = 2L)
    five <- toy_design(stats::setNames(
        lapply(c(12, 8, 10, 3, 6), seq_len), paste0("Path", 1:5)
    ))
    integer <- toy_design(list(Ctrl = c(0, 3, 9), Heat = c(0, 3, 9)), 3L)
    storage.mode(integer) <- "integer"
    list(
        one_path = toy_design(list(A = c(1, 2, 3, 5))),
        two_paths = two_paths,
        replicated = toy_design(list(A = 1:4, B = 1:4), reps = 2L),
        tricky = tricky,
        interleaved = toy_design(list(P1 = 1:5, P2 = 1:4, P3 = 1:6),
            interleave = TRUE
        ),
        five = five,
        integer = integer
    )
}

masigpro_design <- function(e, degree) {
    maSigPro::make.design.matrix(e,
        degree = degree, time.col = 1,
        repl.col = 2, group.cols = 3:ncol(e)
    )
}

test_that(".poly_design() matches maSigPro::make.design.matrix()", {
    skip_if_not_installed("maSigPro")
    e <- new.env()
    utils::data("edesign.abiotic", "NBdesign", package = "maSigPro", envir = e)
    designs <- c(
        oracle_designs(),
        list(abiotic = e$edesign.abiotic, nb = e$NBdesign)
    )
    for (nm in names(designs)) {
        for (degree in 1:5) {
            expect_identical(
                .poly_design(designs[[nm]], degree),
                masigpro_design(designs[[nm]], degree),
                label = paste(nm, "degree", degree)
            )
        }
    }
    # The designs sc.set.poly() builds for the bundled datasets
    for (ob in list(splat_design(1), splat_design(3))) {
        ma <- masigpro_design(pathAssign(ob), ob@Parameters@poly_degree)
        expect_identical(predictors(ob), ma$dis)
        expect_identical(ob@Design@groups.vector, ma$groups.vector)
    }
})

test_that(".path_coeffs() matches maSigPro::reg.coeffs()", {
    skip_if_not_installed("maSigPro")
    for (e in oracle_designs()) {
        gv <- .poly_design(e, degree = 3)$groups.vector
        p <- length(gv)
        coefs <- list(
            seq_len(p + 1) / 7,
            replace(sin(seq_len(p + 1)), c(2, p + 1), NA),
            c(rep(0, p), 1)
        )
        paths <- c(colnames(e)[-(1:2)], "absent")
        for (gvx in list(c(colnames(e)[3], gv), gv)) {
            new <- old <- list()
            for (cf in coefs) {
                for (g in paths) {
                    ma <- maSigPro::reg.coeffs(cf,
                        groups.vector = gvx, group = g
                    )
                    new <- c(new, list(.path_coeffs(cf, gvx, g)))
                    old <- c(old, list(ma))
                }
            }
            expect_identical(new, old)
        }
    }
})
