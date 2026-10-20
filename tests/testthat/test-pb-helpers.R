# Internal binning helpers. extract_interval() counts points in equal-width
# bins and parses the bounds back from the cut() labels, so the bounds are
# rounded to 3 significant digits (ledger S3, kept until Phase 5).
# estBinSize() gives the number of bins of each rule; Freedman-Diaconis and
# Scott give a bin width, turned into a number of bins by the range (ledger g).

test_that("extract_interval() counts equal-width bins over the range", {
    x <- c(0, 0.1, 0.25, 0.5, 0.75, 1)
    bins <- extract_interval(x,
        nBins = 3, bin = "bin", bin.size = "size", lbound = "lo", ubound = "hi"
    )
    expect_identical(colnames(bins), c("lo", "hi", "size"))
    expect_identical(bins$size, c(3, 1, 2))
    # breaks 0, 1/3, 2/3, 1, read back from "[0,0.333]", "(0.333,0.667]", ...
    expect_identical(bins$lo, c(0, 0.333, 0.667))
    expect_identical(bins$hi, c(0.333, 0.667, 1))
    # the first bin is closed on both sides, the others on the right only
    two <- extract_interval(c(2, 4, 4, 9),
        nBins = 2, bin = "bin", bin.size = "size", lbound = "lo", ubound = "hi"
    )
    expect_identical(two$size, c(3, 1))
    expect_identical(two$lo, c(2, 5.5))
})

test_that(".skewness() is the type-3 sample skewness", {
    x <- c(1, 2, 3, 10)
    d <- x - mean(x)
    expect_equal(.skewness(x), mean(d^3) / stats::sd(x)^3)
    expect_identical(.skewness(c(1, 2, 3)), 0)
    expect_identical(.skewness(c(1, NA, 3)), NA)
    expect_true(is.nan(.skewness(c(5, 5, 5))))
})

test_that("estBinSize() counts Freedman-Diaconis and Scott bins (ledger g)", {
    # 0:7: n^(1/3) = 2 and the interquartile range is 3.5, so the
    # Freedman-Diaconis width is 3.5 and the range of 7 holds 2 bins (the
    # width itself was returned as the number of bins before)
    x <- 0:7
    expect_equal(estBinSize(x, 8, 1, "Freedman.Diaconis"), 2)
    expect_equal(estBinSize(x, 8, 1, "Scott.Normal"), 7 / (3.49 * sd(x) / 2))
    # drop_fac scales the number of bins, as for the other rules
    expect_equal(estBinSize(x, 8, 0.5, "Freedman.Diaconis"), 1)
    expect_equal(
        estBinSize(x, 8, 0.5, "Scott.Normal"),
        0.5 * estBinSize(x, 8, 1, "Scott.Normal")
    )
    expect_equal(estBinSize(x, 8, 0.5, "Sturges"), 2)
    # The number of bins does not depend on the units of the pseudotime
    for (method in c("Freedman.Diaconis", "Scott.Normal")) {
        expect_equal(
            estBinSize(x / 7, 8, 1, method), estBinSize(x, 8, 1, method),
            info = method
        )
    }
    # Two points: the interquartile range is half their distance
    expect_equal(estBinSize(c(1, 2), 2, 1, "Freedman.Diaconis"), 2^(1 / 3))
})

test_that("estBinSize() gives one bin when the width is 0 (ledger g)", {
    msg <- function(spread, n, value = "0") {
        paste0(
            "The ", spread, " of the ", n, " pseudotime value(s) is ", value
        )
    }
    iqr <- "interquartile range"
    sdev <- "standard deviation"
    for (method in c("Freedman.Diaconis", "Scott.Normal")) {
        expect_message(
            expect_identical(estBinSize(rep(5, 8), 8, 1, method), 1),
            paste0(
                msg(if (method == "Scott.Normal") sdev else iqr, 8),
                ", so bin_method = '", method,
                "' gives no bin width; they are counted as one bin."
            ),
            fixed = TRUE
        )
    }
    # An interquartile range of 0 with a positive range: Scott has a width
    x <- c(rep(1, 10), 2, 3)
    expect_message(
        expect_identical(estBinSize(x, 12, 1, "Freedman.Diaconis"), 1),
        msg(iqr, 12),
        fixed = TRUE
    )
    expect_silent(estBinSize(x, 12, 1, "Scott.Normal"))
    # One point: an interquartile range of 0, no standard deviation
    expect_message(
        expect_identical(estBinSize(3, 1, 1, "Freedman.Diaconis"), 1),
        msg(iqr, 1),
        fixed = TRUE
    )
    expect_message(
        expect_identical(estBinSize(3, 1, 1, "Scott.Normal"), 1),
        msg(sdev, 1, "undefined"),
        fixed = TRUE
    )
    # drop_fac scales the one bin like any number of bins
    expect_identical(
        suppressMessages(estBinSize(rep(5, 8), 8, 0.5, "Scott.Normal")), 0.5
    )
    # Missing or infinite pseudotime has no number of bins (sc.squeeze()
    # cannot bin it with any rule) and no message
    expect_silent(expect_identical(
        estBinSize(c(1, NA, 3), 3, 1, "Scott.Normal"), NA_real_
    ))
    expect_silent(expect_identical(
        estBinSize(c(1, 2, Inf), 3, 1, "Freedman.Diaconis"), NA_real_
    ))
    expect_error(
        estBinSize(c(1, NA, 3), 3, 1, "Freedman.Diaconis"), "missing values"
    )
})

test_that("sc.squeeze() bins by the Freedman-Diaconis and Scott counts", {
    bins_per_path <- function(ob) {
        as.vector(table(cDense(ob)[[ob@Parameters@path_col]]))
    }
    # splat.sim, pseudotime 1 to 100: both rules give 5 bins per path, where
    # using the width as the number of bins made 21 and 24 (Freedman-Diaconis)
    # and 20 and 23 (Scott); the same on a pseudotime from 0 to 1, where it
    # made one bin per path
    scaled <- splat_scmp()
    cd <- cSparse(scaled)
    cd$Pseudotime <- (cd$Pseudotime - 1) / 99
    cSparse(scaled) <- cd
    for (method in c("Freedman.Diaconis", "Scott.Normal")) {
        ob <- sc.squeeze(splat_scmp(), bin_method = method)
        expect_identical(bins_per_path(ob), c(5L, 5L), info = method)
        expect_identical(ob@Parameters@bin_method, method)
        ob <- sc.squeeze(scaled, bin_method = method)
        expect_identical(bins_per_path(ob), c(5L, 5L), info = method)
    }
    # A path whose pseudotime has an interquartile range of 0 is one bin;
    # the other, 40 evenly spread points (width 1 / 40^(1/3), 3.4 bins) has 4
    cd <- data.frame(
        Pseudotime = c(rep(0.5, 30), seq(0, 1, length.out = 10), (0:39) / 39),
        Path = rep(c("Path1", "Path2"), each = 40),
        row.names = paste0("cell", 1:80)
    )
    counts <- matrix(1, 2, 80, dimnames = list(c("g1", "g2"), rownames(cd)))
    ob <- create_scmp(
        counts = counts, cell_data = cd, ptime_col = "Pseudotime",
        path_col = "Path"
    )
    expect_message(
        ob <- sc.squeeze(ob, bin_method = "Freedman.Diaconis"),
        "The interquartile range of the 40 pseudotime value(s) is 0",
        fixed = TRUE
    )
    expect_identical(bins_per_path(ob), c(1L, 4L))
})
