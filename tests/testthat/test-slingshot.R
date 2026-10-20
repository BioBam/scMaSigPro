# From a slingshot trajectory to the fitted model: slingshot's pseudotimes
# and curve weights go through assign_paths(), create_scmp(), sc.squeeze(),
# sc.restruct() and the fits. The bundled multi.lin.sce holds the t-SNE
# embedding and the Monocle3 branch labels of the same cells.

test_that("slingshot lineages run through assign_paths() and the workflow", {
    skip_if_not_installed("slingshot")
    sce <- load_dataset("multi.lin.sce")
    labels <- SummarizedExperiment::colData(sce)$group
    leaves <- c("path_Y_15", "path_Y_18", "path_Y_52")
    sds <- slingshot::slingshot(
        SingleCellExperiment::reducedDim(sce, "TSNE"),
        clusterLabels = labels, start.clus = "root", end.clus = leaves
    )
    # Name each lineage after the cluster it ends in
    ends <- vapply(slingshot::slingLineages(sds), function(l) {
        l[length(l)]
    }, character(1))
    expect_setequal(ends, leaves)
    time <- slingshot::slingPseudotime(sds)
    weights <- slingshot::slingCurveWeights(sds)
    colnames(time) <- colnames(weights) <- ends
    cells <- assign_paths(time, weights, shared_sep = "|", verbose = FALSE)

    # Every cell is on a lineage, and the labels follow the Monocle3 ones:
    # the root cells are root cells there, and most cells of each end path
    # are on the same branch there
    expect_identical(rownames(cells), colnames(sce))
    expect_gt(mean(labels[cells$Path == "root"] == "root"), 0.95)
    for (leaf in leaves) {
        expect_gt(mean(labels[cells$Path == leaf] == leaf), 0.75)
    }
    link <- paste(intersect(ends, c("path_Y_15", "path_Y_18")),
        collapse = "|"
    )
    expect_true(link %in% cells$Path)

    # The workflow on 40 genes, with the root and the Y_15-Y_18 link
    # restructured into the three end paths
    genes <- rownames(sce)[seq_len(40)]
    counts <- as.matrix(SummarizedExperiment::assay(sce, "counts"))
    ob <- create_scmp(counts[genes, rownames(cells)], cells,
        ptime_col = "Pseudotime", path_col = "Path"
    )
    ob <- suppressMessages(sc.squeeze(ob, verbose = FALSE))
    ob <- sc.restruct(ob,
        end_node_list = as.list(leaves), root_node = "root",
        link_node_list = list(link), link_sep = "|", verbose = FALSE
    )
    expect_setequal(cDense(ob)$Path, leaves)
    ob <- sc.set.poly(ob, poly_degree = 2)
    ob <- sc.p.vector(ob, verbose = FALSE)
    ob <- sc.t.fit(ob, verbose = FALSE)
    expect_valid_scmp(ob)
    expect_gt(length(ob@Profile@non_flat), 0)
    expect_identical(rownames(scmpResults(ob)), genes)
})
