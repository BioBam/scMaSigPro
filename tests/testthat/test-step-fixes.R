# Wrong results and crashes inside the pipeline steps (roadmap #18b), one
# block per ledger item. Exact messages, no snapshot.

# The design of `ob` with its columns renamed
rename_design <- function(ob, names) {
    m <- predictors(ob)
    colnames(m) <- names
    predictors(ob) <- m
    ob
}

test_that("sc.filter(vars = 'each') gives one list per term (ledger F1)", {
    ob <- splat_fitted("default")
    sol <- showSol(ob)
    terms <- sub("p.valor_", "", colnames(sol)[-(1:2)], fixed = TRUE)
    kept <- rownames(sol)[sol[["R-squared"]] >= 0.7 & sol[["p-value"]] <= 0.05]
    # Every term selects the same genes: apply() returned a matrix
    each <- sc.filter(ob, vars = "each", term_p_value = 1)@Significant@genes
    expect_identical(names(each), terms)
    expect_identical(unname(each), rep(list(kept), length(terms)))
    # No term selects a gene
    none <- sc.filter(ob, vars = "each", term_p_value = 0)@Significant@genes
    expect_identical(none, setNames(rep(list(character(0)), 6), terms))
    expect_identical(
        sc.filter(ob, vars = "each", rsq = 1.1)@Significant@genes, none
    )
    # Control: terms that select different numbers of genes
    some <- sc.filter(ob, vars = "each")@Significant@genes
    expect_identical(
        lengths(some), setNames(c(31L, 9L, 14L, 23L, 8L, 13L), terms)
    )
})

test_that("sc.t.fit() uses its p_value argument (ledger b)", {
    pv <- splat_fitted("default")
    strict <- sc.t.fit(pv, p_value = 1e-4, verbose = FALSE)
    expect_identical(strict@Parameters@p_value, 1e-4)
    terms <- as.matrix(showSol(strict))[, -(1:2)]
    expect_true(all(terms <= 1e-4, na.rm = TRUE))
    expect_false(identical(strict@Estimate, pv@Estimate))
    # Control: the stored level is the default
    same <- sc.t.fit(pv, p_value = pv@Parameters@p_value, verbose = FALSE)
    expect_identical(same@Estimate, pv@Estimate)
    # At level 0 no gene keeps a coefficient
    expect_error(sc.t.fit(pv, p_value = 0, verbose = FALSE),
        "No gene has a coefficient significant at p_value = 0.",
        fixed = TRUE
    )
})

test_that("sc.p.vector() records a run without genes (ledger B8)", {
    design <- splat_design(2)
    expect_message(
        none <- sc.p.vector(design, p_value = 1e-300, verbose = FALSE),
        "No significant genes detected. Try changing parameters.",
        fixed = TRUE
    )
    expect_identical(none@Profile@non_flat, character(0))
    expect_identical(none@Profile@p_values, splat_fitted()@Profile@p_values)
    expect_identical(none@Parameters@p_value, 1e-300)
    expect_true(validObject(none))
    # A second run replaces the profile of the first, so sc.t.fit() stops
    rerun <- suppressMessages(
        sc.p.vector(splat_fitted(), p_value = 1e-300, verbose = FALSE)
    )
    expect_identical(rerun@Profile, none@Profile)
    expect_error(sc.t.fit(rerun, verbose = FALSE),
        "No gene with a non-flat profile found",
        fixed = TRUE
    )
})

test_that("sc.p.vector() fits one gene left by the filters (ledger Z3)", {
    design <- splat_design(2)
    counts <- eDense(design)
    gene <- rownames(counts)[1]
    full <- sc.p.vector(design, offset = FALSE, verbose = FALSE)
    ob <- design
    # The other genes have only zero counts, then only missing counts
    for (other in c(0, NA)) {
        x <- counts
        x[-1, ] <- other
        suppressWarnings(eDense(ob, "bulk.counts") <- x)
        one <- suppressMessages(
            sc.p.vector(ob, offset = FALSE, verbose = FALSE)
        )
        expect_identical(one@Profile@p_values, full@Profile@p_values[gene])
    }
})

test_that("the fits do not depend on the design names (ledger N1)", {
    design <- splat_design(2)
    p <- ncol(predictors(design))
    ref <- list(backward = splat_fitted("default"))
    ref$two.ways.forward <- sc.t.fit(ref$backward,
        selection_method = "two.ways.forward", verbose = FALSE
    )
    for (nm in list(
        paste0("term ", seq_len(p), "|x"),
        rep("t", p),
        replace(colnames(predictors(design)), 2, "y"),
        replace(colnames(predictors(design)), 2, "(Intercept)"),
        paste0("t", seq_len(p))
    )) {
        pv <- sc.p.vector(rename_design(design, nm), verbose = FALSE)
        expect_identical(pv@Profile, ref$backward@Profile, info = nm[2])
        for (m in names(ref)) {
            tf <- sc.t.fit(pv, selection_method = m, verbose = FALSE)
            info <- paste(nm[2], m)
            expect_identical(
                unname(showSol(tf)), unname(showSol(ref[[m]])),
                info = info
            )
            expect_identical(
                unname(showCoeff(tf)), unname(showCoeff(ref[[m]])),
                info = info
            )
            expect_identical(colnames(showCoeff(tf)),
                c("beta0", paste0("beta", nm)),
                info = info
            )
        }
    }
})

test_that("Mclust names one-component clusters by gene (ledger C1)", {
    skip_if_not_installed("mclust")
    ob <- splat_filtered()
    genes <- rownames(showCoeff(ob))
    genes <- genes[genes %in% Reduce(union, ob@Significant@genes)]
    clusters <- sc.cluster.trend(ob,
        geneSet = "union", k = 1, cluster_method = "Mclust"
    )@Significant@clusters
    expect_identical(clusters, as.list(setNames(rep(1, length(genes)), genes)))
})

test_that("hclust with another distance clusters the genes (ledger S4)", {
    ob <- splat_filtered()
    genes <- Reduce(union, ob@Significant@genes)
    coeff <- as.matrix(showCoeff(ob))
    coeff <- coeff[rownames(coeff) %in% genes, , drop = FALSE]
    coeff[is.na(coeff)] <- 0
    for (method in c("euclidean", "manhattan")) {
        expected <- cutree(
            hclust(dist(coeff, method = method), method = "ward.D"),
            k = 4
        )
        clusters <- sc.cluster.trend(ob,
            geneSet = "union", k = 4, hclust.distance = method
        )@Significant@clusters
        expect_identical(clusters, as.list(expected), info = method)
    }
    # Control: the default correlation distance
    expect_identical(
        sc.cluster.trend(ob, geneSet = "union", k = 4)@Significant@clusters,
        splat_clustered()@Significant@clusters
    )
})

test_that("split_bins stops when splitting makes no progress (ledger Q2)", {
    # Path1 has 20 of its 60 cells at 0.3: the bin that holds them is too
    # large, and its cells have an interquartile range of 0, so
    # Freedman-Diaconis counts them as one bin and the split gives the same
    # bin back. The loop never ended: the time limit guards a regression
    cd <- data.frame(
        Pseudotime = c(
            seq(0, 1, length.out = 40), rep(0.3, 20), seq(0, 1, length.out = 60)
        ),
        Path = rep(c("Path1", "Path2"), each = 60),
        row.names = paste0("cell", 1:120)
    )
    counts <- matrix(1, 2, 120, dimnames = list(c("g1", "g2"), rownames(cd)))
    ob <- create_scmp(
        counts = counts, cell_data = cd, ptime_col = "Pseudotime",
        path_col = "Path"
    )
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
    setTimeLimit(elapsed = 30, transient = TRUE)
    expect_error(
        suppressMessages(sc.squeeze(ob,
            bin_method = "Freedman.Diaconis", split_bins = TRUE
        )),
        paste(
            "Splitting the bins of 'Path1' with bin_method =",
            "'Freedman.Diaconis' makes no progress"
        ),
        fixed = TRUE
    )
    setTimeLimit(elapsed = Inf)
    # Freedman-Diaconis and Scott once split each large bin of splat.sim into
    # one bin (ledger g); they now make progress
    for (method in c("Freedman.Diaconis", "Scott.Normal")) {
        split <- sc.squeeze(splat_scmp(),
            bin_method = method, split_bins = TRUE
        )
        expect_identical(
            as.vector(table(cDense(split)[[split@Parameters@path_col]])),
            c(7L, 8L),
            info = method
        )
    }
    # Control: a rule whose splits make progress
    split <- sc.squeeze(splat_scmp(), bin_method = "Sturges", split_bins = TRUE)
    expect_identical(ncol(split@Dense), 25L)
})

restruct_multi <- function(ob, link, sep) {
    sc.restruct(ob,
        end_node_list = list("Y_15", "Y_18", "Y_52"), root_node = "Root",
        link_node_list = list(link), link_sep = sep, verbose = FALSE
    )
}

test_that("sc.restruct() uses the path column of the object (ledger i)", {
    x <- multi_cell_data()
    names(x$cell_data)[names(x$cell_data) == "group"] <- "lineage"
    ob <- sc.squeeze(create_scmp(
        counts = x$counts, cell_data = x$cell_data,
        ptime_col = "Monocle3_Pseudotime", path_col = "lineage"
    ))
    rs <- restruct_multi(ob, "Y_18_links_Y_15", "_links_")
    ref <- multi_restruct()
    expect_identical(eDense(rs), eDense(ref))
    bins <- cDense(rs)
    names(bins)[names(bins) == "lineage"] <- "group"
    expect_identical(bins, cDense(ref))
})

test_that("sc.restruct() splits on a fixed link_sep (ledger B5)", {
    x <- multi_cell_data()
    x$cell_data$group <- gsub("_links_", "|", x$cell_data$group, fixed = TRUE)
    ob <- sc.squeeze(create_scmp(
        counts = x$counts, cell_data = x$cell_data,
        ptime_col = "Monocle3_Pseudotime", path_col = "group"
    ))
    fixed <- restruct_multi(ob, "Y_18|Y_15", "|")
    expect_identical(eDense(fixed), eDense(multi_restruct()))
    # Control: the escaped regular expression
    expect_identical(restruct_multi(ob, "Y_18|Y_15", "\\|"), fixed)
})

test_that("sc.t.fit() keeps the coefficients of the genes of sol (ledger T2)", {
    # At this level one gene keeps no significant coefficient
    tf <- sc.t.fit(splat_fitted("default"), p_value = 1e-10, verbose = FALSE)
    genes <- rownames(showSol(tf))
    expect_length(genes, 48L)
    expect_identical(rownames(showCoeff(tf)), genes)
    expect_identical(rownames(showTS(tf)), genes)
    expect_identical(rownames(showGroupCoeff(tf)), genes)
})

test_that("sc.p.vector() and sc.t.fit() fit with the link (ledger S6)", {
    design <- splat_design(2)
    log_fit <- sc.p.vector(design, offset = FALSE, verbose = FALSE)
    # glm() warns for genes whose identity-link mean becomes negative
    id_fit <- suppressWarnings(sc.p.vector(design,
        link = "identity", offset = FALSE, verbose = FALSE
    ))
    family <- id_fit@Parameters@distribution
    expect_identical(id_fit@Parameters@link, "identity")
    expect_identical(family$link, "identity")
    expect_identical(family$family, "Negative Binomial(10)")
    expect_identical(family$linkfun(c(0.5, 2)), c(0.5, 2))
    expect_false(identical(id_fit@Profile@p_values, log_fit@Profile@p_values))
    # The same fits as with the family built by hand
    by_hand <- suppressWarnings(sc.p.vector(design,
        family = MASS::negative.binomial(theta = 10, link = "identity"),
        offset = FALSE, verbose = FALSE
    ))
    expect_identical(by_hand@Profile, id_fit@Profile)
    expect_identical(by_hand@Parameters@link, "identity")
    # sc.t.fit() fits with the family and the link of sc.p.vector()
    id_tf <- suppressWarnings(sc.t.fit(id_fit, verbose = FALSE))
    log_tf <- sc.t.fit(log_fit, verbose = FALSE)
    expect_identical(id_tf@Parameters@link, "identity")
    expect_identical(id_tf@Parameters@distribution$linkfun(2), 2)
    expect_false(identical(showCoeff(id_tf), showCoeff(log_tf)))
    # Control: the default fits have the log link
    default <- splat_fitted("default")
    expect_identical(default@Parameters@link, "log")
    expect_identical(default@Parameters@distribution$linkfun(2), log(2))
    expect_identical(
        sc.p.vector(design, link = "log", verbose = FALSE)@Profile,
        default@Profile
    )
})

test_that("a family given to the fits keeps its own link (ledger S6)", {
    design <- splat_design(2)
    # Recorded as "log", gaussian() needed positive counts to start; no gene
    # is non-flat with it on these counts
    gauss <- suppressMessages(sc.p.vector(design,
        family = gaussian(), offset = FALSE, verbose = FALSE
    ))
    expect_identical(gauss@Parameters@link, "identity")
    expect_identical(gauss@Parameters@distribution$link, "identity")
    expect_identical(
        suppressMessages(sc.p.vector(design,
            family = gaussian(), link = "identity", offset = FALSE,
            verbose = FALSE
        ))@Profile,
        gauss@Profile
    )
    # A family with the log link gives the default fits
    nb <- sc.p.vector(design,
        family = MASS::negative.binomial(theta = 10), link = "log",
        verbose = FALSE
    )
    expect_identical(nb@Profile, splat_fitted("default")@Profile)
    # sc.t.fit() takes another family with its own link
    tf <- sc.t.fit(splat_fitted("default"),
        family = gaussian(), offset = FALSE, verbose = FALSE
    )
    expect_identical(tf@Parameters@link, "identity")
    expect_identical(tf@Parameters@distribution$family, "gaussian")
})
