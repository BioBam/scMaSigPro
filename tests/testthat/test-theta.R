# Estimated theta for the negative binomial family (decision D16):
# sc.p.vector(theta_estimate = "none" | "gene" | "global"). The argument
# checks are in test-input-checks.R, the parallel backends in
# test-parallel.R.

same <- function(x, y) {
    identical(x, y, ignore.environment = TRUE, ignore.bytecode = TRUE)
}

test_that("theta_estimate = 'none' is the default and changes nothing", {
    ob <- splat_design(2)
    expect_true(same(
        sc.p.vector(ob, verbose = FALSE),
        sc.p.vector(ob, verbose = FALSE, theta_estimate = "none")
    ))
    fam <- sc.p.vector(ob, verbose = FALSE)@Parameters@distribution
    expect_null(fam$theta_estimate)
    expect_identical(fam$family, "Negative Binomial(10)")
})

test_that(".theta_gene() is the maximum likelihood theta of glm.nb()", {
    # A gene with an overdispersed profile, its design and the size factors
    ob <- splat_design(2)
    dis <- .fit_frame(ob@Design@predictor_matrix)
    k <- as.matrix(SummarizedExperiment::assay(ob@Dense, "bulk.counts"))
    k <- k[, rownames(dis)] + 1
    sf <- .size_factors(k)
    fam <- MASS::negative.binomial(10)
    i <- which(rownames(k) == "Gene51")
    est <- .theta_gene(i, k, dis, fam, 1e-10, sf, 100)
    y <- as.numeric(k[i, ])
    nb <- MASS::glm.nb(y ~ . + offset(sf),
        data = dis, init.theta = 10,
        control = stats::glm.control(epsilon = 1e-10, maxit = 100)
    )
    expect_equal(unname(est[["theta"]]), nb$theta, tolerance = 1e-4)
    expect_identical(unname(est[c("at_bound", "converged")]), c(0, 1))
    # Counts no more variable than Poisson counts end at the upper bound
    pois <- .theta_gene(
        1, rbind(c(10, 11, 9, 10, 10, 11, 9, 10)),
        data.frame(x1 = rep(0:1, each = 4)), fam, 1e-8, NULL, 100
    )
    expect_identical(unname(pois[["at_bound"]]), 1)
})

test_that("theta_estimate = 'gene' fits every gene with its own theta", {
    ob <- sc.p.vector(splat_design(2), verbose = FALSE, theta_estimate = "gene")
    fam <- ob@Parameters@distribution
    expect_identical(fam$theta_estimate, "gene")
    expect_identical(names(fam$theta_gene), names(ob@Profile@p_values))
    expect_true(all(fam$theta_gene >= 1e-2 & fam$theta_gene <= 1e6))
    expect_identical(
        showParams(ob)$value[showParams(ob)$parameters == "distribution"],
        "Negative Binomial(theta per gene)"
    )
    # The p-value of a gene is the test with its theta
    g <- "Gene51"
    dis <- .fit_frame(ob@Design@predictor_matrix)
    k <- as.matrix(SummarizedExperiment::assay(ob@Dense, "bulk.counts"))
    k <- k[, rownames(dis)] + 1
    p <- .pvector_gene(
        which(rownames(k) == g), k, dis,
        MASS::negative.binomial(fam$theta_gene[[g]]), 1e-8,
        .size_factors(k), 100
    )
    expect_equal(ob@Profile@p_values[[g]], p)
    # sc.t.fit() and plotDiagnostics() use the same theta
    tf <- sc.t.fit(ob, verbose = FALSE)
    expect_identical(tf@Parameters@distribution$theta_gene, fam$theta_gene)
    expect_identical(
        .gene_family(fam, g)$family,
        MASS::negative.binomial(fam$theta_gene[[g]])$family
    )
})

test_that("theta_estimate = 'global' fits every gene with the median", {
    gene <- sc.p.vector(splat_design(2),
        verbose = FALSE, theta_estimate = "gene"
    )
    glob <- sc.p.vector(splat_design(2),
        verbose = FALSE, theta_estimate = "global"
    )
    fam <- glob@Parameters@distribution
    th <- median(gene@Parameters@distribution$theta_gene)
    expect_identical(fam$theta_estimate, "global")
    expect_identical(fam$family, MASS::negative.binomial(th)$family)
    expect_equal(environment(fam$variance)$.Theta, th)
    expect_identical(.gene_family(fam, "Gene51"), fam)
    tf <- sc.t.fit(glob, verbose = FALSE)
    expect_s4_class(tf, "ScMaSigPro")
})

test_that("an estimated theta keeps the link of the family", {
    design <- splat_design(2)
    id <- suppressWarnings(sc.p.vector(design,
        link = "identity", offset = FALSE, theta_estimate = "global",
        verbose = FALSE
    ))
    fam <- id@Parameters@distribution
    expect_identical(fam$link, "identity")
    expect_identical(fam$linkfun(c(0.5, 2)), c(0.5, 2))
    expect_identical(id@Parameters@link, "identity")
})

test_that("theta is estimated on no gene when every gene is left out", {
    design <- splat_design(2)
    none <- suppressMessages(sc.p.vector(design,
        min_counts = Inf, theta_estimate = "global", verbose = FALSE
    ))
    expect_identical(none@Parameters@g, 0L)
    fam <- none@Parameters@distribution
    expect_identical(fam$theta_estimate, "global")
    expect_length(fam$theta_gene, 0L)
    expect_identical(fam$family, "Negative Binomial(10)")
})
