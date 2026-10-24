# The test of each gene's full model in sc.p.vector(): the chi-squared test
# of the deviance difference (test = "chisq", the default) or the F test of
# anova(test = "F"). The argument checks are in test-input-checks.R.

# A hand-made gene: two paths of four bins, a degree-1 design (the path
# dummy, time and their product) and counts that rise on the second path
hand_dis <- function() {
    path <- rep(0:1, each = 4)
    time <- rep(1:4, 2)
    data.frame(x1 = path, x2 = time, x3 = time * path)
}
hand_counts <- rbind(c(16, 12, 21, 14, 15, 22, 27, 41))
hand_offset <- log(c(0.8, 1.1, 0.9, 1.2, 1, 0.9, 1.1, 1))

test_that("test = 'F' is the F test of anova() on the deviance difference", {
    dis <- hand_dis()
    fam <- MASS::negative.binomial(10)
    y <- as.numeric(hand_counts[1, ])
    full <- glm(y ~ .,
        data = dis, family = fam, offset = hand_offset, epsilon = 1e-8,
        maxit = 100
    )
    null <- glm(y ~ 1,
        family = fam, offset = hand_offset, epsilon = 1e-8, maxit = 100
    )
    p_f <- .pvector_gene(1, hand_counts, dis, fam, 1e-8, hand_offset, 100,
        f_test = TRUE
    )
    p_chisq <- .pvector_gene(1, hand_counts, dis, fam, 1e-8, hand_offset, 100)
    expect_identical(p_f, anova(null, full, test = "F")[2, "Pr(>F)"])
    expect_identical(p_chisq, anova(null, full, test = "Chisq")[2, "Pr(>Chi)"])
    # Both divide the deviance difference by the dispersion of the full
    # model; the chi-squared test compares it with the chi-squared
    # distribution, the F test, per degree of freedom, with the F
    # distribution on the residual degrees of freedom of the full model
    df <- null$df.residual - full$df.residual
    stat <- (null$deviance - full$deviance) / summary(full)$dispersion
    expect_identical(c(df, full$df.residual), c(3L, 4L))
    expect_equal(p_chisq, pchisq(stat, df, lower.tail = FALSE),
        tolerance = 1e-12
    )
    expect_equal(p_f, pf(stat / df, df, full$df.residual, lower.tail = FALSE),
        tolerance = 1e-12
    )
    # This gene is significant at 0.05 by the chi-squared test only
    expect_lt(p_chisq, 0.05)
    expect_gt(p_f, 0.05)
})

test_that("a fit without residual degrees of freedom has a p-value of 1", {
    # Four bins and four coefficients
    bins <- c(1, 4, 5, 8)
    dis <- hand_dis()[bins, ]
    counts <- hand_counts[, bins, drop = FALSE]
    fam <- MASS::negative.binomial(10)
    for (f_test in c(FALSE, TRUE)) {
        expect_identical(
            .pvector_gene(1, counts, dis, fam, 1e-8, NULL, 100,
                f_test = f_test
            ),
            1
        )
    }
    # The same with sc.p.vector(): three bins per path and a degree-2 design
    binned <- suppressMessages(sc.squeeze(splat_scmp(), drop_fac = 0.3))
    design <- sc.set.poly(binned, poly_degree = 2)
    expect_identical(dim(predictors(design)), c(6L, 5L))
    fits <- lapply(c(chisq = "chisq", F = "F"), function(test) {
        expect_message(
            ob <- sc.p.vector(design, test = test, verbose = FALSE),
            "No significant genes detected."
        )
        ob@Profile@p_values
    })
    expect_true(all(fits$F == 1))
    expect_identical(fits$F, fits$chisq)
})

test_that("sc.p.vector() gives every gene the p-value of its test", {
    design <- splat_design(2)
    chisq <- sc.p.vector(design, verbose = FALSE)
    f <- sc.p.vector(design, test = "F", verbose = FALSE)
    # The default is the chi-squared test, as before
    expect_identical(
        sc.p.vector(design, test = "chisq", verbose = FALSE)@Profile,
        chisq@Profile
    )
    expect_identical(chisq@Profile, splat_fitted("default")@Profile)
    # The same genes and offsets, each with the p-value of the F test
    expect_identical(names(f@Profile@p_values), names(chisq@Profile@p_values))
    expect_identical(f@Design@offset, chisq@Design@offset)
    dis <- .fit_frame(predictors(design))
    counts <- eDense(design)[, rownames(dis)]
    fam <- MASS::negative.binomial(10)
    for (gene in c("Gene11", "Gene40")) {
        expect_identical(
            f@Profile@p_values[[gene]],
            .pvector_gene(which(rownames(counts) == gene), counts, dis, fam,
                1e-8, f@Design@offset, 100,
                f_test = TRUE
            )
        )
    }
    expect_identical(
        f@Profile@adj_p_values,
        p.adjust(f@Profile@p_values, method = "BH")
    )
})

test_that("the Gaussian family has the F test with either value", {
    dis <- hand_dis()
    p <- vapply(c(FALSE, TRUE), function(f_test) {
        .pvector_gene(1, hand_counts, dis, gaussian(), 1e-8, NULL, 100,
            f_test = f_test
        )
    }, numeric(1))
    y <- as.numeric(hand_counts[1, ])
    full <- glm(y ~ ., data = dis, family = gaussian())
    expect_identical(p, rep(anova(glm(y ~ 1), full, test = "F")[2, 6], 2))
    # On splat.sim no gene is significant after the correction, and every
    # gene has the same p-value with both values
    design <- splat_design(2)
    gauss <- function(test) {
        suppressMessages(sc.p.vector(design,
            family = gaussian(), offset = FALSE, test = test,
            verbose = FALSE
        ))@Profile
    }
    expect_identical(gauss("F"), gauss("chisq"))
})

test_that("the test is kept in the test parameter", {
    design <- splat_design(2)
    expect_identical(design@Parameters@test, "chisq")
    expect_identical(splat_fitted("default")@Parameters@test, "chisq")
    f <- sc.p.vector(design, test = "F", verbose = FALSE)
    expect_identical(f@Parameters@test, "F")
    # sc.t.fit() keeps it, and showParams() lists it
    f <- sc.t.fit(f, verbose = FALSE)
    expect_identical(f@Parameters@test, "F")
    params <- showParams(f)
    expect_identical(params$value[params$parameters == "test"], "F")
})
