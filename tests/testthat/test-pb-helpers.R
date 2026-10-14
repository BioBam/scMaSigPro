# Internal binning helpers. extract_interval() counts points in equal-width
# bins and parses the bounds back from the cut() labels, so the bounds are
# rounded to 3 significant digits (ledger S3, kept until Phase 5).

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
