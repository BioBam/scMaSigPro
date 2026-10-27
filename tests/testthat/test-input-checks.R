# Checks at the start of the exported functions: the object
# they are given, then their arguments. Exact messages, no snapshot.

# Every exported function that takes an ScMaSigPro object, called on `ob`
calls_on <- function(ob, gene) {
    list(
        sc.squeeze = function() sc.squeeze(ob),
        sc.restruct = function() {
            sc.restruct(ob,
                end_node_list = list("Path1"), root_node = "Path2",
                link_node_list = list(), verbose = FALSE
            )
        },
        pb_counts = function() pb_counts(ob),
        sc.set.poly = function() sc.set.poly(ob),
        sc.p.vector = function() sc.p.vector(ob, verbose = FALSE),
        sc.t.fit = function() sc.t.fit(ob, verbose = FALSE),
        sc.filter = function() sc.filter(ob, vars = "groups"),
        sc.cluster.trend = function() sc.cluster.trend(ob),
        showCoeff = function() showCoeff(ob),
        showInflu = function() showInflu(ob),
        showTS = function() showTS(ob),
        showSol = function() showSol(ob),
        showSigProf = function() showSigProf(ob),
        showPoly = function() showPoly(ob),
        showParams = function() showParams(ob),
        showGroupCoeff = function() showGroupCoeff(ob),
        queryCoeff = function() queryCoeff(ob, verbose = FALSE),
        plotBinTile = function() plotBinTile(ob),
        plotIntersect = function() plotIntersect(ob, return = TRUE),
        plotTrend = function() plotTrend(ob, feature_id = gene),
        plotTrendCluster = function() plotTrendCluster(ob, verbose = FALSE),
        plotDiagnostics = function() plotDiagnostics(ob, feature_id = gene)
    )
}

test_that("an object of another class is refused by every function", {
    for (f in calls_on(42, "Gene1")) {
        expect_error(f(), "Please provide object of class 'ScMaSigPro'",
            fixed = TRUE
        )
    }
})

test_that("an object saved before the version slot stops before any work", {
    old <- splat_clustered()
    attr(old, "version") <- NULL
    for (f in calls_on(old, splat_sig_gene())) {
        expect_error(f(), "update it with updateObject() first.",
            fixed = TRUE
        )
    }
})

test_that("a Dense slot without counts asks for pb_counts()", {
    msg <- paste(
        "The 'Dense' slot has no 'bulk.counts' assay: run pb_counts()",
        "after changing the bins with cDense<-."
    )
    ob <- splat_design(2)
    suppressMessages(cDense(ob) <- cDense(ob)[-1, ])
    expect_error(sc.p.vector(ob, verbose = FALSE), msg, fixed = TRUE)
    ob <- splat_clustered()
    suppressMessages(cDense(ob) <- cDense(ob)[-1, ])
    gene <- splat_sig_gene()
    expect_error(sc.t.fit(ob, verbose = FALSE), msg, fixed = TRUE)
    expect_error(plotTrend(ob, feature_id = gene), msg, fixed = TRUE)
    expect_error(showSigProf(ob, includeInflu = TRUE), msg, fixed = TRUE)
    expect_error(plotDiagnostics(ob, feature_id = gene), msg, fixed = TRUE)
})

test_that("the fits need counts for the bins of the design", {
    msg <- "The 'bulk.counts' assay has no counts for some bins of the design"
    # Bins renamed with cDense<- before sc.set.poly(), without pb_counts()
    ob <- splat_binned()
    bins <- cDense(ob)
    rownames(bins) <- paste0("bin", seq_len(nrow(bins)))
    cDense(ob) <- bins
    ob <- sc.set.poly(ob, poly_degree = 2)
    expect_error(sc.p.vector(ob, verbose = FALSE), msg, fixed = TRUE)
    # pb_counts() gives the counts the new names
    expect_identical(
        sc.p.vector(pb_counts(ob), verbose = FALSE)@Profile,
        splat_fitted("default")@Profile
    )
    # Renamed after sc.set.poly(): the design and the counts still agree
    ob <- splat_design(2)
    cDense(ob) <- bins
    expect_identical(
        sc.p.vector(ob, verbose = FALSE)@Profile,
        splat_fitted("default")@Profile
    )
})

test_that("the steps say which earlier step is missing", {
    expect_error(sc.p.vector(splat_binned(), verbose = FALSE),
        "No polynomial design found. Please run 'sc.set.poly()', first.",
        fixed = TRUE
    )
    msg <- "No gene with a non-flat profile found. Please run 'sc.p.vector()'"
    expect_error(sc.t.fit(splat_design(2), verbose = FALSE), msg,
        fixed = TRUE
    )
    # sc.p.vector() that selected no gene
    none <- suppressMessages(
        sc.p.vector(splat_design(2), p_value = 1e-300, verbose = FALSE)
    )
    expect_error(sc.t.fit(none, verbose = FALSE), msg, fixed = TRUE)
    expect_error(sc.cluster.trend(splat_fitted("default")),
        "'Significant' slot is empty, please run 'sc.filter()' first.",
        fixed = TRUE
    )
})

test_that("a choice argument takes one of its values, once", {
    two <- function(x) c(x, x)
    ob <- splat_clustered()
    gene <- splat_sig_gene()
    bins <- "Available binning methods are 'Freedman.Diaconis', 'Sqrt',"
    expect_error(sc.squeeze(splat_scmp(), bin_method = two("Sturges")), bins,
        fixed = TRUE
    )
    agg <- "Invalid 'aggregate' value. Please choose either 'mean' or 'sum'."
    expect_error(sc.squeeze(splat_scmp(), aggregate = "max"), agg,
        fixed = TRUE
    )
    expect_error(pb_counts(splat_binned(), cluster_count_by = two("sum")),
        agg,
        fixed = TRUE
    )
    expect_error(
        sc.restruct(multi_binned(),
            end_node_list = list("Y_15", "Y_18", "Y_52"), root_node = "Root",
            link_node_list = list("Y_18_links_Y_15"), aggregate = "max",
            verbose = FALSE
        ),
        agg,
        fixed = TRUE
    )
    expect_error(
        sc.t.fit(splat_fitted("default"),
            selection_method = c("backward", "forward"), verbose = FALSE
        ),
        "stepwise method must be one of backward, forward,",
        fixed = TRUE
    )
    fitted <- splat_fitted("default")
    expect_error(sc.filter(fitted, vars = c("all", "each")),
        "Invalid selection for 'vars'",
        fixed = TRUE
    )
    expect_error(sc.filter(fitted, vars = "groups", intercept = "some"),
        "Invalid selection for 'intercept'",
        fixed = TRUE
    )
    # intercept is used with vars = "groups" only
    expect_identical(
        sc.filter(fitted, vars = "all", intercept = "some"),
        sc.filter(fitted, vars = "all")
    )
    filtered <- splat_filtered()
    expect_error(
        sc.cluster.trend(filtered, geneSet = c("intersect", "union")),
        "does not exist. Please use one of",
        fixed = TRUE
    )
    for (arg in list(
        list(cluster_by = two("coeff")), list(use_dim = two("col")),
        list(fill_na = two("zero"))
    )) {
        expect_error(
            do.call(sc.cluster.trend, c(list(filtered), arg)),
            "is not a valid option. Please use one of",
            fixed = TRUE
        )
    }
    expect_error(
        sc.cluster.trend(filtered, cluster_method = c("hclust", "kmeans")),
        "is not a valid method. Please use one of hclust, kmeans, Mclust",
        fixed = TRUE
    )
    expect_error(
        plotTrend(ob, feature_id = gene, summary_mode = two("mean")),
        "is not a valid option. Please use one of median, mean",
        fixed = TRUE
    )
    expect_error(plotTrend(ob, feature_id = gene, logType = "ln"),
        "'logType' should be one of 'log2', 'log10', 'log'",
        fixed = TRUE
    )
    # logType is used with logs = TRUE only
    expect_s3_class(
        plotTrend(ob, feature_id = gene, logs = FALSE, logType = "ln"),
        "ggplot"
    )
    expect_error(plotDiagnostics(ob, feature_id = gene, model = two("full")),
        "Invalid 'model', please select one of",
        fixed = TRUE
    )
    expect_error(queryCoeff(ob, query = two("path"), verbose = FALSE),
        "Invalid query, please select one of",
        fixed = TRUE
    )
    expect_error(
        queryCoeff(ob, change = c("increasing", "decreasing"), verbose = FALSE),
        "Invalid change, please select one of",
        fixed = TRUE
    )
})

test_that("the family, link and offset of the fits agree", {
    design <- splat_design(2)
    fitted <- splat_fitted("default")
    # glm() warns for genes whose identity-link mean becomes negative
    identity <- suppressWarnings(sc.p.vector(design,
        link = "identity", offset = FALSE, verbose = FALSE
    ))
    # The size-factor offset needs the log link
    offset_msg <- function(link) {
        paste0(
            "The size-factor offset applies to the log link only: with the '",
            link, "' link, pass offset = FALSE."
        )
    }
    expect_error(sc.p.vector(design, link = "identity", verbose = FALSE),
        offset_msg("identity"),
        fixed = TRUE
    )
    expect_error(sc.p.vector(design, family = gaussian(), verbose = FALSE),
        offset_msg("identity"),
        fixed = TRUE
    )
    expect_error(
        sc.p.vector(design, family = poisson(link = "sqrt"), verbose = FALSE),
        offset_msg("sqrt"),
        fixed = TRUE
    )
    expect_error(sc.t.fit(fitted, family = gaussian(), verbose = FALSE),
        offset_msg("identity"),
        fixed = TRUE
    )
    expect_error(sc.t.fit(identity, offset = TRUE, verbose = FALSE),
        offset_msg("identity"),
        fixed = TRUE
    )
    # A given link must be the link of the family
    expect_error(
        sc.p.vector(design,
            family = MASS::negative.binomial(theta = 10), link = "identity",
            offset = FALSE, verbose = FALSE
        ),
        "'family' has the 'log' link, not link = 'identity'.",
        fixed = TRUE
    )
    expect_error(
        sc.p.vector(design,
            family = gaussian(), link = "log", offset = FALSE, verbose = FALSE
        ),
        "'family' has the 'identity' link, not link = 'log'.",
        fixed = TRUE
    )
    expect_error(sc.t.fit(fitted, link = "identity", verbose = FALSE),
        "'family' has the 'log' link, not link = 'identity'.",
        fixed = TRUE
    )
    expect_error(sc.t.fit(identity, link = "log", verbose = FALSE),
        "'family' has the 'identity' link, not link = 'log'.",
        fixed = TRUE
    )
    # The link is one of two values, given once, and the family an object
    for (bad in list("probit", c("log", "identity"))) {
        msg <- "link function should be either 'log' or 'identity'"
        expect_error(sc.p.vector(design, link = bad, verbose = FALSE), msg,
            fixed = TRUE
        )
        expect_error(sc.t.fit(fitted, link = bad, verbose = FALSE), msg,
            fixed = TRUE
        )
    }
    for (bad in list("gaussian", gaussian)) {
        expect_error(sc.p.vector(design, family = bad, verbose = FALSE),
            "'family' must be a family object",
            fixed = TRUE
        )
    }
    # Control: the stored link given again
    expect_identical(
        sc.t.fit(fitted, link = "log", verbose = FALSE)@Estimate,
        fitted@Estimate
    )
})

test_that("sc.filter() uses 'all' by default", {
    fitted <- splat_fitted("default")
    expect_identical(sc.filter(fitted), sc.filter(fitted, vars = "all"))
    expect_named(sc.filter(fitted)@Significant@genes, "all")
})

test_that("column names and labels are single strings", {
    sce <- load_dataset("splat.sim")
    for (arg in c("path_prefix", "root_label", "ptime_col", "path_col")) {
        for (bad in list(NA_character_, c("a", "b"), 1, character())) {
            args <- stats::setNames(
                list(sce, FALSE, bad), c("object", "verbose", arg)
            )
            expect_error(do.call(as_scmp, args),
                paste0("'", arg, "' must be a single character string."),
                fixed = TRUE
            )
        }
    }
    x <- multi_cell_data()
    expect_error(
        create_scmp(
            counts = x$counts, cell_data = x$cell_data,
            ptime_col = c("Monocle3_Pseudotime", "group"), path_col = "group"
        ),
        "'ptime_col' must be a single character string.",
        fixed = TRUE
    )
    expect_error(sc.squeeze(splat_scmp(), bin_col = c("a", "b")),
        "'bin_col' must be a single character string.",
        fixed = TRUE
    )
    expect_error(sc.squeeze(splat_scmp(), assay_name = c("counts", "counts")),
        "'assay_name' must be a single character string.",
        fixed = TRUE
    )
    expect_error(sc.set.poly(splat_binned(), path_col = NA_character_),
        "'path_col' must be a single character string.",
        fixed = TRUE
    )
    expect_error(pb_counts(splat_binned(), bin_mem_col = 1),
        "'bin_mem_col' must be a single character string.",
        fixed = TRUE
    )
})

test_that("additional_params can be NULL", {
    ob <- splat_scmp()
    expect_identical(sc.squeeze(ob, additional_params = NULL), splat_binned())
    expect_identical(sc.squeeze(ob, additional_params = list()), splat_binned())
    expect_error(sc.squeeze(ob, additional_params = list(TRUE)),
        "Allowed additional parameters are 'use_unique_time_points'.",
        fixed = TRUE
    )
    sce <- load_dataset("splat.sim")
    expect_identical(
        as_scmp(sce, verbose = FALSE, additional_params = NULL),
        as_scmp(sce, verbose = FALSE)
    )
    expect_identical(
        as_scmp(sce, verbose = FALSE, additional_params = list()),
        as_scmp(sce, verbose = FALSE)
    )
})

test_that("sc.cluster.trend() needs enough genes for k", {
    ob <- splat_filtered()
    set <- function(n) {
        genes <- sort(unique(unlist(ob@Significant@genes)))[seq_len(n)]
        ob@Significant <- new("Significant", genes = list(set = genes))
        ob
    }
    msg <- function(n, method, k, need) {
        paste0(
            "The 'set' gene set has ", n, " gene(s), too few for ",
            "cluster_method = '", method, "' with k = ", k,
            ": it needs at least ", need, "."
        )
    }
    expect_error(sc.cluster.trend(set(1), geneSet = "set", k = 1),
        msg(1, "hclust", 1, 2),
        fixed = TRUE
    )
    expect_error(sc.cluster.trend(set(3), geneSet = "set", k = 4),
        msg(3, "hclust", 4, 4),
        fixed = TRUE
    )
    expect_error(
        sc.cluster.trend(set(3),
            geneSet = "set", k = 3,
            cluster_method = "kmeans"
        ),
        msg(3, "kmeans", 3, 4),
        fixed = TRUE
    )
    # Enough genes: as before
    expect_named(
        sc.cluster.trend(set(3), geneSet = "set", k = 3)@Significant@clusters,
        sort(unique(unlist(ob@Significant@genes)))[1:3]
    )
})

test_that("sc.set.poly() takes one of the paths as reference", {
    binned <- splat_binned()
    expect_error(sc.set.poly(binned, reference = "Path3"),
        "'reference' must be one of the paths: Path1, Path2.",
        fixed = TRUE
    )
    for (bad in list(c("Path1", "Path2"), NA_character_, 1, character())) {
        expect_error(sc.set.poly(binned, reference = bad),
            "'reference' must be a single character string.",
            fixed = TRUE
        )
    }
    # Control: the first path given as reference is the default
    expect_identical(
        sc.set.poly(binned, reference = "Path1")@Design,
        splat_design(2)@Design
    )
})

test_that("sc.p.vector() checks the count filter", {
    design <- splat_design(2)
    nbins <- nrow(predictors(design))
    msg <- paste0(
        "'min_bins_expressed' must be a single number from 0 to ", nbins,
        ", the number of bins of the design."
    )
    for (bad in list(nbins + 1, -1, NA, c(1, 2), "1")) {
        expect_error(
            sc.p.vector(design, min_bins_expressed = bad, verbose = FALSE),
            msg,
            fixed = TRUE
        )
    }
    for (bad in list(-1, NA_real_, c(1, 2), "1", TRUE)) {
        expect_error(sc.p.vector(design, min_counts = bad, verbose = FALSE),
            "'min_counts' must be a single number of at least 0.",
            fixed = TRUE
        )
    }
    # Control: the largest threshold is the number of bins
    expect_s4_class(
        sc.p.vector(design, min_bins_expressed = nbins, verbose = FALSE),
        "ScMaSigPro"
    )
})

test_that("an object saved before the count filter stops before any work", {
    old <- splat_clustered()
    params <- old@Parameters
    attr(params, "min_counts") <- NULL
    attr(params, "min_bins_expressed") <- NULL
    old@Parameters <- params
    for (f in calls_on(old, splat_sig_gene())) {
        expect_error(f(), "update it with updateObject() first.",
            fixed = TRUE
        )
    }
    expect_output(print(old), "update it with updateObject()", fixed = TRUE)
})

test_that("test is 'chisq' or 'F', and 'F' needs an estimated dispersion", {
    design <- splat_design(2)
    for (bad in list("f", "Chisq", "LRT", c("chisq", "F"), NA_character_, 1)) {
        expect_error(sc.p.vector(design, test = bad, verbose = FALSE),
            "'test' must be either 'chisq' or 'F'",
            fixed = TRUE
        )
    }
    declared <- MASS::negative.binomial(10)
    declared$dispersion <- 1
    families <- list(poisson(), binomial(link = "log"), declared)
    for (fam in families) {
        expect_error(
            sc.p.vector(design, family = fam, test = "F", verbose = FALSE),
            paste0(
                "test = 'F' needs a family whose dispersion is estimated: ",
                "the '", fam$family, "' family fixes it at 1. ",
                "Use test = 'chisq'."
            ),
            fixed = TRUE
        )
    }
    # Control: the chi-squared test takes the Poisson family
    expect_s4_class(
        sc.p.vector(design, family = poisson(), verbose = FALSE),
        "ScMaSigPro"
    )
})

test_that("an object saved before the test parameter stops before any work", {
    old <- splat_clustered()
    params <- old@Parameters
    attr(params, "test") <- NULL
    old@Parameters <- params
    for (f in calls_on(old, splat_sig_gene())) {
        expect_error(f(), "update it with updateObject() first.",
            fixed = TRUE
        )
    }
    expect_output(print(old), "update it with updateObject()", fixed = TRUE)
})

test_that("theta_estimate is one of three values", {
    design <- splat_design(2)
    for (bad in list("each", c("gene", "global"), NA_character_, 1)) {
        expect_error(
            sc.p.vector(design, theta_estimate = bad, verbose = FALSE),
            "'theta_estimate' must be one of 'none', 'gene' or 'global'",
            fixed = TRUE
        )
    }
    expect_error(
        sc.p.vector(design,
            theta_estimate = "gene", family = gaussian(), offset = FALSE,
            verbose = FALSE
        ),
        "'theta_estimate' needs a negative binomial family",
        fixed = TRUE
    )
})

test_that("a relabelled family is refused", {
    design <- splat_design(2)
    msg <- "'family' is labelled with the 'identity' link but has the link"
    nb <- MASS::negative.binomial(10)
    nb$link <- "identity"
    expect_error(
        sc.p.vector(design, family = nb, offset = FALSE, verbose = FALSE),
        msg,
        fixed = TRUE
    )
    expect_error(sc.t.fit(splat_relabelled(), verbose = FALSE), msg,
        fixed = TRUE
    )
    expect_error(
        sc.t.fit(splat_fitted("no_offset"), family = nb, verbose = FALSE),
        msg,
        fixed = TRUE
    )
    # The other direction: identity functions under the log label
    gauss <- gaussian()
    gauss$link <- "log"
    expect_error(
        sc.p.vector(design, family = gauss, offset = FALSE, verbose = FALSE),
        "'family' is labelled with the 'log' link but has the link functions",
        fixed = TRUE
    )
    # Controls: families with their own link, and a link make.link() does
    # not know, which is not checked
    expect_identical(
        .true_link(MASS::negative.binomial(10, link = "sqrt")), "sqrt"
    )
    expect_identical(.true_link(poisson(link = power(1 / 3))), NA_character_)
    expect_s4_class(
        sc.p.vector(design,
            family = poisson(link = power(1 / 3)), offset = FALSE,
            verbose = FALSE
        ),
        "ScMaSigPro"
    )
})

test_that("sc.squeeze() takes any drop_fac of at least 0.3", {
    # Srivastava et al. (2024), Section 2.1: the drop factor is at least
    # 0.3, with no upper bound; the message said "between 0.3 and 1"
    expect_error(sc.squeeze(splat_scmp(), drop_fac = 0.29),
        "Invalid value for 'drop_fac'. It should be at least 0.3.",
        fixed = TRUE
    )
    # Sturges gives the 97 and 103 cells of the two paths 7.6 and 7.7 bins
    for (k in c(0.3, 1.5)) {
        ob <- sc.squeeze(splat_scmp(), drop_fac = k)
        expect_valid_scmp(ob)
        expect_identical(
            as.vector(table(cDense(ob)[[ob@Parameters@path_col]])),
            if (k == 0.3) c(3L, 3L) else c(12L, 12L),
            info = k
        )
    }
})

test_that("arguments without effect are gone", {
    gone <- list(
        plotTrendCluster = c("plot", "includeInflu"),
        plotIntersect = c("width_ratio", "show_sets_size"),
        queryCoeff = "includeInflu",
        sc.cluster.trend = "verbose"
    )
    for (fn in names(gone)) {
        expect_false(any(gone[[fn]] %in% names(formals(fn))), label = fn)
    }
})
