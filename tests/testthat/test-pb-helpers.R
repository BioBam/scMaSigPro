# Internal binning helpers. extract_interval() counts points in equal-width
# bins whose bounds are the numeric breaks (ledger S3), and .assign_bins()
# puts points in bins with the same intervals.
# estBinSize() gives the number of bins of each rule; Freedman-Diaconis and
# Scott give a bin width, turned into a number of bins by the range (ledger g).

test_that("extract_interval() counts equal-width bins over the range", {
    x <- c(0, 0.1, 0.25, 0.5, 0.75, 1)
    bins <- extract_interval(x,
        nBins = 3, bin.size = "size", lbound = "lo", ubound = "hi"
    )
    expect_identical(colnames(bins), c("lo", "hi", "size"))
    expect_identical(bins$size, c(3, 1, 2))
    # the bounds are the breaks 0, 1/3, 2/3 and 1, not the 0.333 and 0.667 of
    # their cut() labels
    breaks <- seq(0, 1, length.out = 4)
    expect_identical(bins$lo, breaks[-4])
    expect_identical(bins$hi, breaks[-1])
    two <- extract_interval(c(2, 4, 4, 9),
        nBins = 2, bin.size = "size", lbound = "lo", ubound = "hi"
    )
    expect_identical(two$size, c(3, 1))
    expect_identical(two$lo, c(2, 5.5))
    # the last bin is closed on both sides, the others on the left only: the
    # point on the break 2 is counted in the second bin
    on_break <- extract_interval(0:4,
        nBins = 2, bin.size = "size", lbound = "lo", ubound = "hi"
    )
    expect_identical(on_break$size, c(2, 3))
})

test_that(".assign_bins() puts each point in its bin", {
    # [0, 2), [2, 4) and [4, 6], the last closed: a point on the bound of two
    # bins goes to the upper one, as with cut(right = FALSE,
    # include.lowest = TRUE)
    x <- 0:6
    bin <- .assign_bins(x, c(0, 2, 4), c(2, 4, 6))
    expect_identical(bin, c(1L, 1L, 2L, 2L, 3L, 3L, 3L))
    expect_identical(bin, cut(x, c(0, 2, 4, 6),
        right = FALSE, include.lowest = TRUE, labels = FALSE
    ))
    # bins with a gap between them, and points outside every bin
    expect_identical(
        .assign_bins(c(-1, 1, 2.5, 3, 7), c(0, 3), c(2, 6)),
        c(NA, 1L, NA, 2L, NA)
    )
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

test_that("estBinSize() rounds the Sturges bins up before drop_fac", {
    # Srivastava et al. (2024), Section 2.1: the number of bins is the
    # Sturges formula, [log2(N)+1], times the drop factor k. 200 points:
    # log2(200) + 1 = 8.64 is rounded up to 9 bins, and drop_fac 1.5 gives
    # 13.5, which sc.squeeze() rounds up to 14; it rounded 1.5 x 8.64 =
    # 12.97 up to 13
    x <- seq_len(200)
    expect_identical(estBinSize(x, 200, 1, "Sturges"), 9)
    expect_identical(estBinSize(x, 200, 1.5, "Sturges"), 13.5)
    # The other rules are not rounded before drop_fac
    expect_identical(estBinSize(x, 200, 1.5, "Sqrt"), 1.5 * sqrt(200))
    expect_identical(estBinSize(x, 200, 1.5, "Rice"), 1.5 * 2 * 200^(1 / 3))
    # multi.lin.sce, Sturges, drop_fac 1.5: the 676 cells of Y_18 give
    # log2(676) + 1 = 10.40 bins, rounded up to 11, times 1.5: 16.5, so 17
    # bins of equal width (16 before)
    ob <- sc.squeeze(multi_scmp(), drop_fac = 1.5)
    bins <- cDense(ob)
    y18 <- bins[bins[[ob@Parameters@path_col]] == "Y_18", ]
    width <- (max(y18$scmp_u_bound) - min(y18$scmp_l_bound)) / 17
    expect_equal(
        y18$scmp_u_bound - y18$scmp_l_bound, rep(width, nrow(y18))
    )
    # drop_fac = 1 keeps the number of bins
    expect_identical(splat_binned(), sc.squeeze(splat_scmp(), drop_fac = 1))
})
