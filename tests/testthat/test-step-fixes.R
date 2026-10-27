# Wrong results and crashes inside the pipeline steps, one block per
# fixed defect. Exact messages, no snapshot.

# The design of `ob` with its columns renamed
rename_design <- function(ob, names) {
    m <- predictors(ob)
    colnames(m) <- names
    predictors(ob) <- m
    ob
}

test_that("sc.filter(vars = 'each') gives one list per term", {
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
        lengths(some), setNames(c(39L, 9L, 17L, 32L, 12L, 10L), terms)
    )
})

test_that("sc.t.fit() uses its p_value argument", {
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

test_that("sc.p.vector() records a run without genes", {
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

test_that("sc.p.vector() fits one gene left by the filters", {
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

test_that("the fits do not depend on the design names", {
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

test_that("Mclust names one-component clusters by gene", {
    skip_if_not_installed("mclust")
    ob <- splat_filtered()
    genes <- rownames(showCoeff(ob))
    genes <- genes[genes %in% Reduce(union, ob@Significant@genes)]
    clusters <- sc.cluster.trend(ob,
        geneSet = "union", k = 1, cluster_method = "Mclust"
    )@Significant@clusters
    expect_identical(clusters, as.list(setNames(rep(1, length(genes)), genes)))
})

test_that("hclust with another distance clusters the genes", {
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

test_that("split_bins stops when splitting makes no progress", {
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
    # one bin; they now make progress
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
    expect_identical(ncol(split@Dense), 29L)
})

restruct_multi <- function(ob, link, sep) {
    sc.restruct(ob,
        end_node_list = list("Y_15", "Y_18", "Y_52"), root_node = "Root",
        link_node_list = list(link), link_sep = sep, verbose = FALSE
    )
}

test_that("sc.restruct() uses the path column of the object", {
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

test_that("sc.restruct() splits on a fixed link_sep", {
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

test_that("sc.t.fit() keeps the coefficients of the genes of sol", {
    # At this level one gene keeps no significant coefficient
    tf <- sc.t.fit(splat_fitted("default"), p_value = 1e-14, verbose = FALSE)
    genes <- rownames(showSol(tf))
    expect_length(genes, 51L)
    expect_identical(rownames(showCoeff(tf)), genes)
    expect_identical(rownames(showTS(tf)), genes)
    expect_identical(rownames(showGroupCoeff(tf)), genes)
})

test_that("sc.p.vector() and sc.t.fit() fit with the link", {
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

test_that("a family given to the fits keeps its own link", {
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

# Three binned paths whose labels sort differently in the C locale ("B",
# "a", "c") and in a case-insensitive collation ("a", "B", "c"); `levels`
# makes the path column a factor with these levels
three_paths_binned <- function(labels = c("a", "B", "c"), levels = NULL) {
    n <- 20
    cells <- paste0("cell", seq_len(n * length(labels)))
    path <- rep(labels, each = n)
    if (!is.null(levels)) {
        path <- factor(path, levels = levels)
    }
    cell_data <- data.frame(
        Pseudotime = rep(seq_len(n), length(labels)),
        Path = path,
        row.names = cells
    )
    # Deterministic counts: gene i rises along path k with slope i * k / 400
    time <- cell_data$Pseudotime
    k <- match(as.character(cell_data$Path), labels)
    counts <- t(vapply(seq_len(10), function(i) {
        round(exp(2 + i * k * time / 400 + sin(i * time) / 4))
    }, numeric(length(cells))))
    dimnames(counts) <- list(paste0("g", 1:10), cells)
    ob <- create_scmp(
        counts = counts, cell_data = cell_data,
        ptime_col = "Pseudotime", path_col = "Path"
    )
    sc.squeeze(ob)
}

# The reference path of a design: its shortest groups.vector label
reference_of <- function(ob) {
    gv <- ob@Design@groups.vector
    gv[nchar(gv) == min(nchar(gv))][1]
}

# `code` evaluated with LC_COLLATE = `locale`; NULL when it is not available
with_collate <- function(locale, code) {
    old <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
    if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) {
        return(NULL)
    }
    force(code)
}

test_that("sc.squeeze() orders the paths in C-locale order", {
    in_c <- with_collate("C", sc.set.poly(three_paths_binned(), 1))
    expect_identical(reference_of(in_c), "B")
    expect_identical(unique(cDense(in_c)$Path), c("B", "a", "c"))
    expect_identical(colnames(predictors(in_c))[1:2], c("avsB", "cvsB"))
    # The same design in a locale that collates case-insensitively
    in_en <- with_collate("en_US.UTF-8", {
        if (identical(order(c("B", "a")), c(2L, 1L))) {
            sc.set.poly(three_paths_binned(), 1)
        }
    })
    skip_if(is.null(in_en), "no case-insensitive en_US.UTF-8 collation")
    expect_identical(in_en@Design, in_c@Design)
    expect_identical(cDense(in_en), cDense(in_c))
})

test_that("a factor path column keeps its level order", {
    binned <- three_paths_binned(levels = c("c", "a", "B"))
    expect_identical(unique(cDense(binned)$Path), c("c", "a", "B"))
    expect_identical(reference_of(sc.set.poly(binned, 1)), "c")
})

test_that("sc.set.poly(reference = ) names the reference path", {
    binned <- three_paths_binned()
    ob <- sc.set.poly(binned, 1, reference = "c")
    expect_identical(reference_of(ob), "c")
    expect_identical(
        colnames(predictors(ob)),
        c(
            "Bvsc", "avsc", "scmp_binned_pseudotime",
            "scmp_binned_pseudotimexB", "scmp_binned_pseudotimexa"
        )
    )
    expect_identical(colnames(pathAssign(ob))[3:5], c("c", "B", "a"))
    # The bins keep their order; only the parameterisation changes
    expect_identical(
        rownames(predictors(ob)), rownames(predictors(sc.set.poly(binned, 1)))
    )
    # The first path given again is the default
    expect_identical(
        sc.set.poly(binned, 1, reference = "B")@Design,
        sc.set.poly(binned, 1)@Design
    )
    # The path polynomials of sc.t.fit() use the same reference
    design <- sc.set.poly(binned, 2, reference = "c")
    tf <- sc.t.fit(sc.p.vector(design, verbose = FALSE), verbose = FALSE)
    expect_identical(tf@Estimate@path[1], "c")
    expect_identical(
        colnames(showGroupCoeff(tf))[1:3], c("c_beta0", "c_beta1", "c_beta2")
    )
})

test_that("full-model p-values do not depend on the reference", {
    binned <- three_paths_binned()
    p_of <- function(ref) {
        ob <- sc.set.poly(binned, 2, reference = ref)
        sc.p.vector(ob, verbose = FALSE)@Profile@p_values
    }
    p_b <- p_of("B")
    expect_equal(p_of("a"), p_b, tolerance = 1e-8)
    expect_equal(p_of("c"), p_b, tolerance = 1e-8)
})

test_that("the forward methods name the influential bins", {
    ob <- splat_fitted("default")
    for (m in c("forward", "two.ways.forward")) {
        inf <- showInflu(sc.t.fit(ob, selection_method = m, verbose = FALSE))
        expect_gt(ncol(inf), 0)
        expect_identical(rownames(inf), rownames(predictors(ob)))
        expect_false(anyNA(inf), info = m)
        expect_true(all(colSums(inf == "dffit/cook.d") > 0), info = m)
    }
})

test_that("the count filter keeps every gene by default", {
    design <- splat_design(2)
    fitted <- splat_fitted("default")
    expect_identical(fitted@Parameters@g, nrow(eDense(design)))
    expect_identical(fitted@Parameters@min_counts, 0)
    expect_identical(fitted@Parameters@min_bins_expressed, 0)
    expect_identical(
        sc.p.vector(design,
            min_counts = 0, min_bins_expressed = 0, verbose = FALSE
        )@Profile,
        fitted@Profile
    )
})

test_that("the count filter leaves out the genes below it", {
    design <- splat_design(2)
    fitted <- sc.p.vector(design, verbose = FALSE)
    counts <- eDense(design)[, rownames(predictors(design))]
    total <- rowSums(counts)
    cut <- sort(total)[11]
    msgs <- capture_messages(
        by_total <- sc.p.vector(design, min_counts = cut)
    )
    expect_match(msgs,
        paste("Leaving out", sum(total < cut), "of 100 genes below"),
        fixed = TRUE, all = FALSE
    )
    kept <- rownames(counts)[total >= cut]
    expect_identical(names(by_total@Profile@p_values), kept)
    expect_identical(by_total@Parameters@g, length(kept))
    expect_identical(by_total@Parameters@min_counts, as.numeric(cut))
    # The size factors are those of every gene, so the kept genes keep their
    # p-values; only the multiple-testing correction changes
    expect_identical(by_total@Design@offset, fitted@Design@offset)
    expect_identical(
        by_total@Profile@p_values, fitted@Profile@p_values[kept]
    )
    expect_identical(
        by_total@Profile@adj_p_values,
        p.adjust(fitted@Profile@p_values[kept], method = "BH")
    )
    # Bins with a count above 0
    expressed <- rowSums(counts > 0)
    nbins <- ncol(counts)
    by_bins <- sc.p.vector(design,
        min_bins_expressed = nbins, verbose = FALSE
    )
    expect_identical(
        names(by_bins@Profile@p_values), rownames(counts)[expressed >= nbins]
    )
    expect_identical(by_bins@Parameters@min_bins_expressed, as.numeric(nbins))
})

test_that("sc.p.vector() fits an object with one gene", {
    design <- splat_design(2)
    gene <- rownames(eDense(design))[1]
    one <- design
    one@Sparse <- one@Sparse[gene, ]
    one@Dense <- one@Dense[gene, ]
    expect_valid_scmp(one)
    full <- sc.p.vector(design, offset = FALSE, verbose = FALSE)
    fit <- sc.p.vector(one, offset = FALSE, verbose = FALSE)
    expect_identical(fit@Profile@p_values, full@Profile@p_values[gene])
    expect_identical(fit@Parameters@g, 1L)
})

test_that("sc.t.fit() refuses a family relabelled by 0.0.5", {
    old <- splat_relabelled()
    expect_error(sc.t.fit(old, verbose = FALSE),
        "run sc.p.vector() on it again.",
        fixed = TRUE
    )
    # The functions are those of the log link, with which the fits were made
    expect_identical(.true_link(old@Parameters@distribution), "log")
    # Running sc.p.vector() again fits with the identity link
    refit <- suppressWarnings(sc.p.vector(old,
        link = "identity", offset = FALSE, verbose = FALSE
    ))
    expect_identical(.true_link(refit@Parameters@distribution), "identity")
    expect_s4_class(
        suppressWarnings(sc.t.fit(refit, verbose = FALSE)), "ScMaSigPro"
    )
})

# The fitted means of the selected model of `gene`, refitted on the columns
# of the design whose term the model kept, one per bin
selected_means <- function(ob, gene) {
    sol <- showSol(ob)[gene, -(1:3)]
    terms <- !is.na(unlist(sol))
    x <- predictors(ob)[, terms, drop = FALSE]
    fit <- stats::glm.fit(cbind(1, x), eDense(ob)[gene, rownames(x)],
        family = ob@Parameters@distribution,
        control = stats::glm.control(
            epsilon = ob@Parameters@epsilon, maxit = ob@Parameters@max_it
        )
    )
    fit$fitted.values
}

test_that("identity-link fits are drawn on the count scale", {
    design <- splat_design(2)
    # glm() warns for genes whose identity-link mean becomes negative
    pv <- suppressWarnings(sc.p.vector(design,
        link = "identity", offset = FALSE, verbose = FALSE
    ))
    id <- suppressWarnings(sc.t.fit(pv, verbose = FALSE))
    gene <- rownames(showSol(id))[1]
    p <- plotTrend(id, gene, significant = FALSE)
    built <- ggplot2::ggplot_build(p)$data
    expect_identical(p$labels$y, "Pseudobulk Expression")
    # The points are the counts, not their logs
    bins <- cDense(id)
    counts <- eDense(id)[gene, rownames(bins)]
    expect_identical(built[[1]]$y, unname(counts))
    # The curve of each path passes through the fitted means of its bins
    means <- selected_means(id, gene)
    curve <- p$layers$curves$data
    for (path in unique(bins$Path)) {
        on_path <- bins$Path == path
        at_bins <- stats::approx(
            curve$x[curve$path == path], curve$y[curve$path == path],
            xout = bins$scmp_binned_pseudotime[on_path]
        )$y
        expect_equal(at_bins, unname(means[rownames(bins)[on_path]]),
            tolerance = 1e-6
        )
    }
    # Explicit logs = TRUE logs the points as before
    logged <- plotTrend(id, gene, significant = FALSE, logs = TRUE)
    expect_identical(
        ggplot2::ggplot_build(logged)$data[[1]]$y, unname(log(counts + 1))
    )
    # A size factor left in the design by an earlier log-link fit is not
    # used for an identity-link fit
    stale <- id
    stale@Design@offset <- splat_fitted("default")@Design@offset
    expect_identical(
        ggplot2::ggplot_build(plotTrend(stale, gene, significant = FALSE))$data,
        built
    )
})

test_that("log-link fits are drawn on the log scale", {
    cl <- splat_clustered()
    gene <- splat_sig_gene()
    default <- plotTrend(cl, gene, lines = TRUE)
    logged <- plotTrend(cl, gene, lines = TRUE, logs = TRUE)
    expect_identical(
        ggplot2::ggplot_build(default)$data,
        ggplot2::ggplot_build(logged)$data
    )
    expect_identical(default$labels$y, "log(Pseudobulk Expression)")
    expect_identical(
        plotTrendCluster(cl, verbose = FALSE)$labels$y,
        "log(Pseudobulk Expression)"
    )
    # Explicit logs = FALSE keeps the label of plotTrendCluster()
    expect_identical(
        plotTrendCluster(cl, logs = FALSE, verbose = FALSE)$labels$y,
        "log(Pseudobulk Expression)"
    )
})

test_that("plotTrendCluster() resolves logs as plotTrend()", {
    design <- splat_design(2)
    pv <- suppressWarnings(sc.p.vector(design,
        link = "identity", offset = FALSE, verbose = FALSE
    ))
    id <- suppressWarnings(sc.t.fit(pv, verbose = FALSE))
    id <- sc.filter(id, rsq = 0)
    id <- sc.cluster.trend(id, geneSet = "all", k = 2)
    p <- plotTrendCluster(id, verbose = FALSE)
    expect_identical(p$labels$y, "Pseudobulk Expression")
    raw <- ggplot2::ggplot_build(p)$data
    expect_identical(
        ggplot2::ggplot_build(plotTrendCluster(id,
            logs = FALSE, verbose = FALSE
        ))$data,
        raw
    )
    expect_identical(
        plotTrendCluster(id, logs = TRUE, verbose = FALSE)$labels$y,
        "log(Pseudobulk Expression)"
    )
})

# sc.squeeze() on two paths with the same pseudotime values p: the bins of
# Path1 in order, with their stored bounds and the pseudotime of their cells
toy_bins <- function(p, ...) {
    cd <- data.frame(
        Pseudotime = c(p, p), Path = rep(c("Path1", "Path2"), each = length(p)),
        row.names = paste0("cell", seq_along(c(p, p)))
    )
    counts <- matrix(1, 2, nrow(cd),
        dimnames = list(c("g1", "g2"), rownames(cd))
    )
    ob <- create_scmp(
        counts = counts, cell_data = cd, ptime_col = "Pseudotime",
        path_col = "Path"
    )
    ob <- suppressMessages(sc.squeeze(ob, ...))
    bins <- cDense(ob)
    bins <- bins[bins[[ob@Parameters@path_col]] == "Path1", ]
    members <- strsplit(bins[[ob@Parameters@bin_mem_col]], "|", fixed = TRUE)
    list(
        ob = ob, lower = bins$scmp_l_bound, upper = bins$scmp_u_bound,
        members = lapply(members, function(m) cd[m, "Pseudotime"])
    )
}

test_that("sc.squeeze() bins the cells by the numeric breaks", {
    # 4 points, Sturges: 3 bins, breaks 0, 1/3, 2/3 and 1. The bounds were
    # read back from the cut() labels "[0,0.333]", "(0.333,0.667]" and
    # "(0.667,1]", which put 0.3332 and 0.6668 in the middle bin
    bins <- toy_bins(c(0, 0.3332, 0.6668, 1))
    expect_identical(bins$members, list(c(0, 0.3332), c(0.6668, 1)))
    expect_equal(bins$lower, c(0, 2 / 3))
    expect_equal(bins$upper, c(1 / 3, 1))
    # multi.lin.sce, Sturges: five cells of Y_52 lie between the label 30.5
    # and the break 30.53546 between its bins 4 and 5
    ob <- multi_binned()
    bins <- cDense(ob)
    members <- strsplit(bins[[ob@Parameters@bin_mem_col]], "|", fixed = TRUE)
    bin_of <- setNames(rep(rownames(bins), lengths(members)), unlist(members))
    cells <- c("C520", "C910", "C966", "C1500", "C1723")
    expect_identical(unname(bin_of[cells]), rep("Y_52_bin_4", 5))
    expect_equal(bins["Y_52_bin_4", "scmp_u_bound"], 30.53546, tolerance = 1e-6)
    # The bins do not depend on the origin of the pseudotime: shifted by
    # 1000, 50 of the 200 cells of splat.sim were in another bin
    shifted <- splat_scmp()
    ptime_col <- shifted@Parameters@ptime_col
    cd <- cSparse(shifted)
    cd[[ptime_col]] <- cd[[ptime_col]] + 1000
    cSparse(shifted) <- cd
    mem_col <- shifted@Parameters@bin_mem_col
    expect_identical(
        cDense(sc.squeeze(shifted))[[mem_col]],
        cDense(splat_binned())[[mem_col]]
    )
})

test_that("split_bins keeps every cell", {
    # 13 points, Sturges: 5 bins of width 2. The bin from 2 to 4 holds 10
    # cells, more than the 6.82 allowed, and is split over the range of its
    # cells; the first new bin read "[2.24,2.58]" from its cut() label, so
    # the cell at 2.2366 was in no bin and was removed from the object
    p <- c(0, 0.5, 2.2366, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 3.95, 10)
    expect_no_warning(bins <- toy_bins(p, split_bins = TRUE))
    expect_identical(ncol(bins$ob@Sparse), 26L)
    expect_identical(bins$members, list(
        c(0, 0.5), c(2.2366, 2.5), c(2.7, 2.9), 3.1, c(3.3, 3.5),
        c(3.7, 3.9, 3.95), 10
    ))
    # The 5 new bins have equal width over 2.2366 to 3.95, the range of the
    # cells, but the first starts at 2 and the last ends at 4: they cover
    # the bin from 2 to 4 without gaps
    expect_identical(bins$lower[c(1, 2, 7)], c(0, 2, 8))
    expect_identical(bins$upper[c(1, 6, 7)], c(2, 4, 10))
    expect_identical(bins$lower[3:6], bins$upper[2:5])
    expect_equal(bins$upper[2:5], seq(2.2366, 3.95, length.out = 6)[2:5])
})

test_that("the bins are closed on the left and open on the right", {
    # Srivastava et al. (2024), Section 2.1: "bins are structured as
    # left-closed and right-open intervals, except for the last bin, which is
    # closed on both ends". A test meant to apply this rule compared a
    # logical vector with max(), so every bin was closed on both ends and a
    # cell on the bound of two bins went to the earlier one.
    # 0 to 10, Sturges: 5 bins, breaks 0, 2, 4, 6, 8 and 10
    bins <- toy_bins(0:10)
    expect_identical(bins$members, list(0:1, 2:3, 4:5, 6:7, 8:10))
    expect_identical(bins$lower, c(0, 2, 4, 6, 8))
    # 7 points, 4 bins, breaks 0, 2.5, 5, 7.5 and 10: 5 is in the third bin
    bins <- toy_bins(c(0, 2, 4, 5, 6, 8, 10))
    expect_identical(bins$members, list(c(0, 2), 4, c(5, 6), c(8, 10)))
    # splat.sim, Sturges: Cell169, at the break 50 of Path1, opens bin 5
    ob <- splat_binned()
    bins <- cDense(ob)
    members <- strsplit(bins[[ob@Parameters@bin_mem_col]], "|", fixed = TRUE)
    names(members) <- rownames(bins)
    expect_true("Cell169" %in% members[["Path1_bin_5"]])
    expect_identical(bins["Path1_bin_5", "scmp_l_bound"], 50)
    # The bins of a split follow the rule too: every cell of splat.sim, whose
    # pseudotime takes whole values, lies in [lower, upper) of its bin or at
    # the end of the last bin of its path
    ob <- sc.squeeze(splat_scmp(), split_bins = TRUE)
    p <- ob@Parameters
    bins <- cDense(ob)
    members <- strsplit(bins[[p@bin_mem_col]], "|", fixed = TRUE)
    bin_of <- rep(seq_len(nrow(bins)), lengths(members))
    ptime <- cSparse(ob)[unlist(members), p@ptime_col]
    lower <- bins$scmp_l_bound[bin_of]
    upper <- bins$scmp_u_bound[bin_of]
    path_end <- ave(bins$scmp_u_bound, bins[[p@path_col]], FUN = max)
    at_end <- ptime == upper & upper == path_end[bin_of]
    expect_true(all(ptime >= lower & (ptime < upper | at_end)))
})

test_that("split_bins keeps a path with a single bin", {
    # Path1 has 30 of its 40 cells at 0.5, an interquartile range of 0, so
    # Freedman-Diaconis gives it one bin. The sizes of one bin have no
    # standard deviation, and the split failed with "missing value where
    # TRUE/FALSE needed"; Path2 has 4 bins of 10 cells
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
    method <- "Freedman.Diaconis"
    split <- suppressMessages(sc.squeeze(ob,
        bin_method = method, split_bins = TRUE
    ))
    nosplit <- suppressMessages(sc.squeeze(ob, bin_method = method))
    expect_identical(cDense(split), cDense(nosplit))
    expect_identical(ncol(split@Sparse), 80L)
})

test_that("split_bins divides the bins above the mean + sd of their path", {
    # Srivastava et al. (2024), Section 2.2: the mean and the standard
    # deviation "of the sizes of all bins per path" establish the maximum,
    # mean + sd, and "Bins larger than the maximum size are divided".
    # splat.sim, Sturges: the bins of Path2 hold 20, 9, 14, 12, 10, 11, 17
    # and 10 cells, mean 12.875 and sd 3.87, so the maximum is 16.75 and the
    # bins of 20 and 17 cells are divided. The mean and the sd were rounded
    # to 13 and 4 before they were added, which allowed 17 cells
    ob <- splat_binned()
    p <- ob@Parameters
    bins <- cDense(ob)
    path2 <- bins[bins[[p@path_col]] == "Path2", ]
    sizes <- path2[[p@bin_size_col]]
    expect_identical(sizes, c(20, 9, 14, 12, 10, 11, 17, 10))
    expect_equal(mean(sizes) + sd(sizes), 16.74568, tolerance = 1e-6)
    expect_identical(path2$scmp_l_bound[7], 75.25)
    expect_identical(path2$scmp_u_bound[7], 87.625)
    split <- sc.squeeze(splat_scmp(), split_bins = TRUE)
    bins <- cDense(split)
    split2 <- bins[bins[[p@path_col]] == "Path2", ]
    expect_identical(
        split2[[p@bin_size_col]],
        c(3, 7, 2, 2, 6, 9, 14, 12, 10, 11, 2, 6, 2, 3, 4, 10)
    )
    # The bin of 17 cells, from 75.25 to 87.625, is now five bins
    in_bin7 <- split2$scmp_l_bound >= 75.25 & split2$scmp_u_bound <= 87.625
    expect_identical(split2[in_bin7, p@bin_size_col], c(2, 6, 2, 3, 4))
    expect_true(all(split2[[p@bin_size_col]] < 16.75))
})

test_that("prune_bins removes the bins below the mean - sd of their path", {
    # Srivastava et al. (2024), Section 2.2: the mean and the standard
    # deviation "of the sizes of all bins per path" establish the minimum,
    # mean - sd, and the bins "smaller than the minimum size can be
    # optionally removed". splat.sim, Sturges: the bins of Path1 hold 14, 15,
    # 7, 17, 13, 11, 9 and 11 cells, minimum 12.125 - 3.27 = 8.85, so the bin
    # of 7 cells goes and the bin of 9 stays; Path2 (minimum 9.004) loses its
    # bin of 9. The minimum was taken over the bins of both paths, 12.5 -
    # 3.48 = 9.02, which also removed the bin of 9 cells of Path1
    ob <- sc.squeeze(splat_scmp(), prune_bins = TRUE)
    p <- ob@Parameters
    bins <- cDense(ob)
    expect_identical(
        split(bins[[p@bin_size_col]], bins[[p@path_col]]),
        list(
            Path1 = c(14, 15, 17, 13, 11, 9, 11),
            Path2 = c(20, 14, 12, 10, 11, 17, 10)
        )
    )
    expect_identical(
        rownames(bins)[bins[[p@path_col]] == "Path1"],
        paste0("Path1_bin_", c(1, 2, 4, 5, 6, 7, 8))
    )
    expect_valid_scmp(ob)
    # multi.lin.sce, Sturges: the bins of Root hold 184, 75, 17, 46, 49, 32,
    # 29, 0, 7 and 14 cells, mean 45.3 and sd 53.6, so the minimum is
    # negative and no bin is removed, and no bin of the other paths is below
    # the minimum of its path. The absolute value of the minimum over all
    # paths, 10.66, removed 11 of the 44 bins
    expect_identical(
        sc.squeeze(multi_scmp(), prune_bins = TRUE), multi_binned()
    )
})

test_that("split_bins and prune_bins set their limits from the first bins", {
    # The limits are set once per path, from the bins of the rule, before
    # any bin is divided. 9 points, Sturges: 5 bins of width 1.26 with 1, 4,
    # 2, 0 and 2 cells, mean 1.8 and sd 1.48, so the maximum is 3.28 and the
    # minimum 0.32. The bin of 4 cells is divided into bins of 2, 0 and 2
    # cells, and the bin of 1 cell is kept: the bins after the division, of
    # 1, 2, 2, 2 and 2 cells, would give a minimum of 1.35 and remove it, as
    # the minimum over the bins of both paths did
    p <- c(0.7, 2.2, 2.3, 2.9, 3.1, 3.4, 3.6, 6.7, 7)
    divided <- list(0.7, c(2.2, 2.3), c(2.9, 3.1), c(3.4, 3.6), c(6.7, 7))
    expect_identical(toy_bins(p, split_bins = TRUE)$members, divided)
    expect_identical(
        toy_bins(p, split_bins = TRUE, prune_bins = TRUE)$members, divided
    )
    # Without the division, the same minimum keeps every bin
    expect_identical(
        toy_bins(p, prune_bins = TRUE)$members,
        list(0.7, c(2.2, 2.3, 2.9, 3.1), c(3.4, 3.6), c(6.7, 7))
    )
})

test_that("the offset is the logarithm of the size factors of the bins", {
    # Srivastava et al. (2024), Section 2.3: the offset "is calculated as the
    # logarithm of size factors for each bin". The size factors are the
    # medians over the genes of the ratios of the raw counts of a bin to
    # their geometric means (Supplementary Section S1.3), over the genes
    # without a zero count. On splat.sim they range from 0.59 to 1.68, so
    # the size factor itself, the offset before, was no normalisation
    design <- splat_design(2)
    fitted <- sc.p.vector(design, verbose = FALSE)
    counts <- eDense(design)[, rownames(predictors(design))]
    log_k <- log(counts[rowSums(counts == 0) == 0, ])
    log_sf <- apply(log_k - rowMeans(log_k), 2, stats::median)
    expect_equal(fitted@Design@offset, log_sf, tolerance = 1e-12)
    expect_identical(fitted@Design@offset, log(.size_factors(counts)))
    expect_equal(
        unname(round(exp(fitted@Design@offset[c(3, 9)]), 4)),
        c(0.5905, 1.6798)
    )
    # The fit of a gene is the glm() of its raw counts with that offset
    # (the non-flat gene with the largest p-value: expect_equal() compares
    # p-values below its tolerance as absolute differences)
    p_non_flat <- fitted@Profile@p_values[fitted@Profile@non_flat]
    gene <- names(which.max(p_non_flat))
    y <- counts[gene, ]
    nb <- MASS::negative.binomial(10)
    full <- glm(y ~ .,
        data = .fit_frame(predictors(design)), family = nb,
        offset = log_sf, epsilon = 1e-8, maxit = 100
    )
    null <- glm(y ~ 1,
        family = nb, offset = log_sf, epsilon = 1e-8, maxit = 100
    )
    expect_equal(
        fitted@Profile@p_values[[gene]],
        anova(null, full, test = "Chisq")[2, 5],
        tolerance = 1e-8
    )
    # log_offset, which chose between the size factor and its logarithm, is
    # no longer an argument
    expect_error(sc.p.vector(design, log_offset = TRUE, verbose = FALSE),
        "unused argument (log_offset = TRUE)",
        fixed = TRUE
    )
    expect_error(sc.t.fit(fitted, log_offset = TRUE, verbose = FALSE),
        "unused argument (log_offset = TRUE)",
        fixed = TRUE
    )
    expect_false(methods::.hasSlot(fitted@Parameters, "log_offset"))
})

test_that("a fit without offset stores none, and the plots use none", {
    # sc.set.poly() stores no offset, and neither does a fit without offset,
    # also after a fit with one, whose offsets stayed in the object
    expect_identical(splat_design(2)@Design@offset, numeric(0))
    pv <- sc.p.vector(splat_fitted("default"), offset = FALSE, verbose = FALSE)
    expect_identical(pv@Design@offset, numeric(0))
    expect_false(pv@Parameters@offset)
    tf <- sc.t.fit(splat_fitted("default"), offset = FALSE, verbose = FALSE)
    expect_identical(tf@Design@offset, numeric(0))
    expect_valid_scmp(tf)
    # So refitting without offset gives the object of a fit without offset,
    # and its plots, whose points are the logs of the counts + 1, not of the
    # counts divided by the size factors of the earlier fit
    refit <- sc.t.fit(pv, verbose = FALSE)
    fresh <- splat_fitted("no_offset")
    expect_identical(refit, fresh)
    gene <- rownames(showSol(fresh))[1]
    built <- ggplot2::ggplot_build(plotTrend(refit, gene, significant = FALSE))
    counts <- eDense(fresh)[gene, rownames(cDense(fresh))]
    expect_identical(built$data[[1]]$y, unname(log(counts + 1)))
    # With the offset, the points are divided by the size factors
    fitted <- splat_fitted("default")
    gene <- rownames(showSol(fitted))[1]
    built <- ggplot2::ggplot_build(plotTrend(fitted, gene, significant = FALSE))
    counts <- eDense(fitted)[gene, rownames(cDense(fitted))]
    expect_equal(
        built$data[[1]]$y,
        unname(log(counts / exp(fitted@Design@offset) + 1)),
        tolerance = 1e-12
    )
})

test_that("sc.t.fit() reuses the offset of sc.p.vector()", {
    # The offset is computed once, by sc.p.vector() from every gene it
    # fitted; sc.t.fit() computed it again from the non-flat genes
    pv <- sc.p.vector(splat_design(2), verbose = FALSE)
    tf <- sc.t.fit(pv, verbose = FALSE)
    expect_identical(tf@Design@offset, pv@Design@offset)
    # The selected model of a gene is the glm() of its raw counts on the
    # selected terms with that offset
    sol <- showSol(tf)
    gene <- rownames(sol)[sol[["p.valor_beta0"]] <= 0.05][1]
    terms <- !is.na(unlist(sol[gene, -(1:3)]))
    x <- predictors(tf)[, terms, drop = FALSE]
    fit <- stats::glm.fit(cbind(1, x), eDense(tf)[gene, rownames(x)],
        family = MASS::negative.binomial(10), offset = pv@Design@offset,
        control = stats::glm.control(epsilon = 1e-8, maxit = 100)
    )
    expect_equal(
        unname(unlist(showCoeff(tf)[gene, c(TRUE, terms)])),
        unname(fit$coefficients),
        tolerance = 1e-8
    )
    # Without a stored offset, sc.t.fit() computes the same one
    none <- sc.p.vector(splat_design(2), offset = FALSE, verbose = FALSE)
    none@Profile <- pv@Profile
    computed <- sc.t.fit(none, offset = TRUE, verbose = FALSE)
    expect_identical(computed@Design@offset, pv@Design@offset)
    expect_identical(computed@Estimate, tf@Estimate)
})

test_that("the size factors need a gene without a zero count", {
    # Srivastava et al. (2024), Supplementary Section S1.3: the size factor
    # of a bin is the median of the ratios of its counts to the geometric
    # means of the genes; a gene with a zero count has a geometric mean of 0
    # and is left out. With one gene left, the size factors are its counts
    # over their geometric mean
    design <- splat_design(2)
    bins <- rownames(predictors(design))
    k <- eDense(design)
    keep <- rownames(k)[rowSums(k[, bins] == 0) == 0][1]
    k[rownames(k) != keep, bins[1]] <- 0
    expect_warning(eDense(design) <- k, "Overwriting assay", fixed = TRUE)
    fitted <- sc.p.vector(design, verbose = FALSE)
    y <- log(k[keep, bins])
    expect_equal(fitted@Design@offset, y - mean(y), tolerance = 1e-12)
    # With no such gene, the fit stops; the count filter does not help, as
    # the size factors are computed before it
    k[keep, bins[1]] <- 0
    expect_warning(eDense(design) <- k, "Overwriting assay", fixed = TRUE)
    msg <- "Cannot compute the size factors of the offset: every gene has"
    expect_error(sc.p.vector(design, verbose = FALSE), msg, fixed = TRUE)
    expect_error(
        sc.p.vector(design, min_bins_expressed = length(bins), verbose = FALSE),
        msg,
        fixed = TRUE
    )
    expect_s4_class(
        sc.p.vector(design, offset = FALSE, verbose = FALSE), "ScMaSigPro"
    )
})
