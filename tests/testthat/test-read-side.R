# The read side (roadmap #18c): accessors, show helpers and plots return what
# they document. One block per ledger item; exact messages, no snapshot.

test_that("eSparse() and eDense() return the assay they name (ledger k)", {
    sce <- load_dataset("splat.sim")
    counts <- SingleCellExperiment::counts(sce)
    SummarizedExperiment::assays(sce) <- list(
        logcounts = log1p(as.matrix(counts)), counts = counts
    )
    ob <- as_scmp(sce,
        verbose = FALSE,
        additional_params = list(
            labels_exist = TRUE, exist_ptime_col = "Step",
            exist_path_col = "Group"
        )
    )
    expect_identical(eSparse(ob), as.matrix(counts))
    ob <- splat_binned()
    bulk <- eDense(ob)
    SummarizedExperiment::assays(ob@Dense) <- list(
        scaled = bulk * 2, bulk.counts = bulk
    )
    expect_identical(eDense(ob), bulk)
    # Control: without a counts assay the default is the first assay
    SummarizedExperiment::assays(sce) <- list(x = counts)
    ob <- as_scmp(sce,
        verbose = FALSE,
        additional_params = list(
            labels_exist = TRUE, exist_ptime_col = "Step",
            exist_path_col = "Group"
        )
    )
    expect_identical(eSparse(ob), as.matrix(counts))
})

test_that("the assay setters take any matrix and a default name (ledger A1)", {
    ob <- splat_binned()
    m <- eSparse(ob)
    expect_warning(eSparse(ob) <- m * 2, "Overwriting assay 'counts'.",
        fixed = TRUE
    )
    expect_identical(eSparse(ob), m * 2)
    expect_message(eSparse(ob, "sparse") <- as(m, "dgCMatrix"),
        "Creating a new assay as 'sparse'.",
        fixed = TRUE
    )
    # The integer counts come back as doubles from the dgCMatrix
    expect_identical(eSparse(ob, "sparse"), m + 0)
    suppressMessages(eSparse(ob, "frame") <- as.data.frame(m))
    expect_identical(eSparse(ob, "frame"), m + 0)
    d <- eDense(ob)
    expect_warning(eDense(ob) <- d + 1, "Overwriting assay 'bulk.counts'.",
        fixed = TRUE
    )
    expect_identical(eDense(ob), d + 1)
    suppressMessages(eDense(ob, "frame") <- as.data.frame(d))
    expect_identical(eDense(ob, "frame"), d)
})

test_that("showInflu() is empty when no gene is influential (ledger S7)", {
    expect_error(showInflu(splat_design(2)),
        "'influential' is not computed yet",
        fixed = TRUE
    )
    # Without the one influential gene of this fit, no gene is influential
    ob <- splat_fitted("no_offset")
    ob@Profile@non_flat <- setdiff(
        ob@Profile@non_flat, colnames(showInflu(ob))
    )
    ob <- suppressMessages(sc.t.fit(ob, verbose = FALSE))
    expect_identical(dim(showInflu(ob)), c(0L, 0L))
    expect_identical(showCoeff(ob, includeInflu = FALSE), showCoeff(ob))
    expect_identical(showSol(ob, includeInflu = FALSE), showSol(ob))
})

test_that("plotDiagnostics() loads patchwork itself (ledger c)", {
    # skip_if_not_installed() would load patchwork, which is the point
    skip_if(!nzchar(system.file(package = "patchwork")), "needs patchwork")
    if (!"package:patchwork" %in% search()) {
        try(unloadNamespace("patchwork"), silent = TRUE)
    }
    p <- plotDiagnostics(splat_filtered(), feature_id = splat_sig_gene())
    expect_s3_class(p, "patchwork")
})

test_that("printing a restructured object lists its bin paths (ledger S8)", {
    out <- capture.output(print(multi_restruct()))
    expect_true("Branching Paths: Y_15, Y_52, Y_18" %in% out)
    # Control: before binning, the paths of the cells
    out <- capture.output(print(splat_scmp()))
    expect_true("Branching Paths: Path2, Path1" %in% out)
})

test_that("plotTrend() uses every power of the degree (ledger d)", {
    ob <- sc.squeeze(splat_scmp(), bin_method = "Sqrt")
    ob <- sc.p.vector(sc.set.poly(ob, poly_degree = 7), verbose = FALSE)
    ob <- sc.filter(sc.t.fit(ob, verbose = FALSE), vars = "groups")
    gene <- sort(unique(unlist(ob@Significant@genes)))[1]
    curve <- plotTrend(ob, gene)$layers$curves$data
    curve <- curve[curve$path == "Path1", ]
    beta <- unlist(showGroupCoeff(ob)[gene, paste0("Path1_beta", 0:7)])
    beta[is.na(beta)] <- 0
    expected <- vapply(curve$x, function(x) sum(beta * x^(0:7)), numeric(1))
    expect_equal(curve$y, expected, tolerance = 1e-10)
})

test_that("plotTrendCluster() names the genes it leaves out (ledger m)", {
    ob <- splat_clustered()
    left_out <- names(ob@Significant@clusters)[1]
    ob@Significant@genes <- lapply(ob@Significant@genes, setdiff, left_out)
    expect_warning(
        plotTrendCluster(ob, significant = TRUE, verbose = FALSE),
        paste("No trend for 1 gene(s), left out of the plot:", left_out),
        fixed = TRUE
    )
    # Control: no warning when every gene has a trend
    expect_no_warning(plotTrendCluster(splat_clustered(), verbose = FALSE))
})

test_that("plotTrendCluster() draws lines with or without offset (ledger L2)", {
    ob <- splat_clustered()
    # With an offset the lines are the points, also without the points layer
    p <- plotTrendCluster(ob, lines = TRUE, verbose = FALSE)
    lines_only <- plotTrendCluster(ob,
        lines = TRUE, points = FALSE, verbose = FALSE
    )
    expect_identical(lines_only$layers[[1]]$data, p$layers[[1]]$data)
    # Without an offset the lines summarise the counts of each path
    ob@Design@offset[] <- 0
    p <- plotTrendCluster(ob, lines = TRUE, verbose = FALSE)
    lines <- p$layers[[2]]$data
    expect_type(lines$y_axis, "double")
    expect_setequal(lines$group, c("Path1", "Path2"))
})

test_that("plotTrendCluster() titles each cluster by its genes (ledger S5)", {
    set.seed(1)
    ob <- sc.cluster.trend(splat_filtered(),
        geneSet = "union", k = 4, cluster_method = "kmeans"
    )
    clusters <- unlist(ob@Significant@clusters)
    p <- plotTrendCluster(ob, verbose = FALSE)
    points <- p$layers[[1]]$data
    for (id in unique(clusters)) {
        genes <- names(clusters)[clusters == id]
        title <- paste0("Cluster: ", id, " (", length(genes), " Features)")
        # The median over the genes of the cluster, per bin and path
        each <- do.call(rbind, lapply(genes, function(g) {
            plotTrend(ob, g, significant = FALSE)$layers$points$data
        }))
        expected <- .summarise_by(
            each, c("pooled.time", "path"), "pb.counts", median
        )
        got <- points[points$cluster == title, ]
        got <- got[order(got$x_axis, got$group), ]
        expect_identical(got$y_axis, expected$pb.counts, info = title)
    }
})
