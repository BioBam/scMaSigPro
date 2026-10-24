# .size_factors() computes median-of-ratios size factors (Anders and Huber
# 2010). The expected values below are worked out by hand from the
# definition: geometric mean of each gene, ratio of every count to it, median
# of the ratios per sample (on the log scale).

counts_4x3 <- function() {
    matrix(c(
        2, 4, 8,
        1, 1, 8,
        3, 9, 1,
        0, 6, 6
    ), nrow = 4, byrow = TRUE, dimnames = list(
        c("g1", "g2", "g3", "g4"), c("s1", "s2", "s3")
    ))
}

test_that(".size_factors() is the median of ratios to the geometric means", {
    # Geometric means: g1 64^(1/3) = 4, g2 8^(1/3) = 2, g3 27^(1/3) = 3; g4 has
    # a zero, so its geometric mean is 0 and it is skipped.
    # Ratios: s1 = 0.5, 0.5, 1; s2 = 1, 0.5, 3; s3 = 2, 4, 1/3.
    expect_equal(.size_factors(counts_4x3()), c(s1 = 0.5, s2 = 1, s3 = 2))
    expect_error(.size_factors(counts_4x3(), "median"), "should be one of")
    # a vector (sc.t.fit() with no non-flat gene) is rejected as before
    expect_error(.size_factors(c(1, 2, 3)), "at least two dimensions")
})

test_that("type = 'poscounts' uses the positive counts of every gene", {
    # g4 = (0, 6, 6): geometric mean = cube root of 6 * 6 (the m-th root over
    # all 3 samples); its ratio 6 / 36^(1/3) enters the medians of s2 and s3
    # only. With four ratios, the median on the log scale is the geometric
    # mean of the middle two.
    r4 <- 6 / 36^(1 / 3)
    expect_equal(
        .size_factors(counts_4x3(), type = "poscounts"),
        c(s1 = 0.5, s2 = sqrt(1 * r4), s3 = sqrt(r4 * 2))
    )
    # Without zeros the two types agree
    k <- counts_4x3()[1:3, ]
    expect_identical(.size_factors(k, "poscounts"), .size_factors(k))
})

test_that("a dgeMatrix (dgCMatrix + 1, as in sc.t.fit) is handled", {
    dense <- Matrix::Matrix(counts_4x3(), sparse = TRUE) + 1
    expect_s4_class(dense, "dgeMatrix")
    expect_identical(.size_factors(dense), .size_factors(counts_4x3() + 1))
})

test_that("genes with only zero counts are skipped", {
    k <- counts_4x3()
    with_zero_gene <- rbind(k, g5 = 0)
    expect_identical(.size_factors(with_zero_gene), .size_factors(k))
    expect_identical(
        .size_factors(with_zero_gene, "poscounts"),
        .size_factors(k, "poscounts")
    )
})

test_that("it stops when no gene can be used; missing counts are skipped", {
    every_gene_has_a_zero <- rbind(c(0, 1, 2), c(3, 0, 1))
    expect_error(
        .size_factors(every_gene_has_a_zero),
        "every gene contains at least one zero or missing count",
        fixed = TRUE
    )
    expect_true(all(.size_factors(every_gene_has_a_zero, "poscounts") > 0))
    expect_error(
        .size_factors(matrix(0, 2, 3), "poscounts"),
        "every gene has only zero counts or a missing count",
        fixed = TRUE
    )
    # A gene with a missing count is skipped, like a gene with a zero
    k <- counts_4x3()[1:3, ]
    with_na <- rbind(k, g4 = c(5, NA, 5))
    expect_identical(.size_factors(with_na), .size_factors(k))
    expect_identical(
        .size_factors(with_na, "poscounts"),
        .size_factors(k, "poscounts")
    )
    # With poscounts, a sample without positive counts gets NA
    sf <- .size_factors(cbind(counts_4x3(), s4 = 0), "poscounts")
    expect_true(is.na(sf[["s4"]]))
    expect_true(all(is.finite(sf[c("s1", "s2", "s3")])))
})

test_that(".size_factors() gives the example of Supplementary Section S1.3", {
    # Srivastava et al. (2024), Supplementary Section S1.3: three genes in
    # three bins. Geometric means: 144.22, 246.62 and 407.16; the size factor
    # of a bin is the median of the ratios of its counts to them, here the
    # ratio of Gene 3 in every bin
    k <- matrix(c(100, 150, 200, 200, 300, 250, 300, 450, 500),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste("Gene", 1:3), paste("Bin", 1:3))
    )
    expect_equal(
        .size_factors(k), k["Gene 3", ] / (300 * 450 * 500)^(1 / 3),
        tolerance = 1e-12
    )
    expect_equal(unname(round(.size_factors(k), 4)), c(0.7368, 1.1052, 1.2280))
})
