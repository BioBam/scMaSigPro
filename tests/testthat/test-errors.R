# Error paths: every case must error, and the message is snapshotted so that
# wording changes are reviewed deliberately.

bad_cells <- data.frame(x = 1:2)

cases <- list(
    "as_scmp rejects an unknown 'from'" =
        quote(as_scmp(load_dataset("splat.sim"), from = "seurat")),
    "as_scmp no longer converts Monocle3 objects (from = 'cds')" =
        quote(as_scmp(load_dataset("splat.sim"), from = "cds")),
    "as_scmp rejects an object that is not a SingleCellExperiment" =
        quote(as_scmp(SummarizedExperiment::SummarizedExperiment(list(counts = matrix(1, 2, 2))))),
    "sc.squeeze rejects a non-ScMaSigPro input" =
        quote(sc.squeeze(42)),
    "sc.squeeze rejects an unknown pseudotime column" =
        quote(sc.squeeze(splat_scmp(), ptime_col = "nope")),
    "sc.squeeze rejects an unknown path column" =
        quote(sc.squeeze(splat_scmp(), path_col = "nope")),
    "sc.squeeze rejects an unknown bin method" =
        quote(sc.squeeze(splat_scmp(), bin_method = "Knuth")),
    "sc.squeeze rejects drop_fac below 0.3" =
        quote(sc.squeeze(splat_scmp(), drop_fac = 0.1)),
    "sc.set.poly rejects a non-ScMaSigPro input" =
        quote(sc.set.poly(list())),
    "sc.p.vector rejects an unknown link" =
        quote(sc.p.vector(splat_design(2), link = "probit", verbose = FALSE)),
    "sc.p.vector refuses the offset with the identity link" =
        quote(sc.p.vector(splat_design(2), link = "identity", verbose = FALSE)),
    "sc.p.vector rejects a link that is not the link of family" =
        quote(sc.p.vector(splat_design(2),
            family = MASS::negative.binomial(theta = 10), link = "identity",
            verbose = FALSE
        )),
    "sc.p.vector rejects a family that is not a family object" =
        quote(sc.p.vector(splat_design(2),
            family = "gaussian", verbose = FALSE
        )),
    "sc.p.vector rejects min_na larger than the number of bins" =
        quote(sc.p.vector(splat_design(2), min_na = 1000, verbose = FALSE)),
    "sc.t.fit rejects an unknown selection method" =
        quote(sc.t.fit(splat_fitted("default"), selection_method = "sideways", verbose = FALSE)),
    "sc.filter rejects an unknown vars" =
        quote(sc.filter(splat_fitted("default"), vars = "nothing")),
    "showSol errors before sc.t.fit" =
        quote(showSol(splat_binned())),
    "showCoeff errors before sc.t.fit" =
        quote(showCoeff(splat_binned())),
    "eSparse rejects an unknown assay" =
        quote(eSparse(splat_binned(), i = "nope")),
    "eDense rejects an unknown assay" =
        quote(eDense(splat_binned(), i = "nope")),
    "plotTrend rejects an unknown feature" =
        quote(plotTrend(splat_clustered(), feature_id = "NotAGene")),
    "plotTrend refuses a non-significant feature by default" =
        quote(plotTrend(splat_clustered(), feature_id = splat_nonsig_gene())),
    "plotTrend needs at least one layer" =
        quote(plotTrend(splat_clustered(), feature_id = splat_sig_gene(), curves = FALSE, lines = FALSE, points = FALSE)),
    "plotBinTile errors before sc.squeeze" =
        quote(plotBinTile(splat_scmp())),
    "plotIntersect errors before sc.filter" =
        quote(plotIntersect(splat_fitted("default"))),
    "plotTrendCluster errors before sc.cluster.trend" =
        quote(plotTrendCluster(splat_filtered())),
    "sc.cluster.trend rejects an unknown gene set" =
        quote(sc.cluster.trend(splat_filtered(), geneSet = "nope")),
    "sc.cluster.trend rejects an unknown method" =
        quote(sc.cluster.trend(splat_filtered(), cluster_method = "dbscan")),
    "queryCoeff rejects an unknown query" =
        quote(queryCoeff(splat_fitted("default"), query = "nope", verbose = FALSE)),
    "queryCoeff rejects an unknown change" =
        quote(queryCoeff(splat_fitted("default"), change = "sideways", verbose = FALSE)),
    "create_scmp rejects mismatched cells" =
        quote(create_scmp(counts = matrix(1, 2, 3), cell_data = bad_cells, ptime_col = "x", path_col = "x")),
    "sc.restruct errors before sc.squeeze" =
        quote(sc.restruct(splat_scmp(), end_node_list = list("Path1"), root_node = "Path1", link_node_list = list()))
)

for (nm in names(cases)) {
    test_that(nm, {
        # Snapshots are skipped when NOT_CRAN is unset (R CMD check on CI and on
        # Bioconductor); the error itself is asserted everywhere.
        expect_error(eval(cases[[nm]]))
        expect_snapshot(error = TRUE, eval(cases[[nm]]))
    })
}
