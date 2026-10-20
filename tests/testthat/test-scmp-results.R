# scmpResults() and scmpToRowData(): the per-gene results of each step in
# one DataFrame, and that table written into the rowData() of an SCE.

test_that("scmpResults() has the columns of the steps that have run", {
    genes <- rownames(load_dataset("splat.sim"))
    res <- scmpResults(splat_scmp())
    expect_s4_class(res, "DataFrame")
    expect_identical(rownames(res), genes)
    expect_identical(ncol(res), 0L)
    pv <- sc.p.vector(splat_design(2), verbose = FALSE)
    expect_identical(
        colnames(scmpResults(pv)), c("p_value", "adj_p_value", "non_flat")
    )
    fitted <- c(
        "p_value", "adj_p_value", "non_flat", "model_p_value", "r_squared",
        "influential"
    )
    expect_identical(colnames(scmpResults(splat_fitted())), fitted)
    lists <- paste0("sig_", names(splat_filtered()@Significant@genes))
    expect_identical(colnames(scmpResults(splat_filtered())), c(fitted, lists))
    expect_identical(
        colnames(scmpResults(splat_clustered())), c(fitted, lists, "cluster")
    )
})

test_that("scmpResults() holds the values of the slots", {
    # The splat.sim workflow with Gene1 and Gene2 set to zero counts, so
    # that sc.p.vector() leaves them out
    sce <- load_dataset("splat.sim")
    counts <- SummarizedExperiment::assay(sce, "counts")
    counts[c("Gene1", "Gene2"), ] <- 0
    SummarizedExperiment::assay(sce, "counts") <- counts
    ob <- as_scmp(sce,
        verbose = FALSE,
        additional_params = list(
            labels_exist = TRUE, exist_ptime_col = "Step",
            exist_path_col = "Group"
        )
    )
    ob <- sc.set.poly(sc.squeeze(ob), poly_degree = 2)
    ob <- sc.t.fit(sc.p.vector(ob, verbose = FALSE), verbose = FALSE)
    ob <- sc.filter(ob, rsq = 0.7, vars = "groups")
    ob <- sc.cluster.trend(ob, geneSet = "union", k = 4)
    res <- scmpResults(ob)
    genes <- rownames(res)
    expect_identical(genes, rownames(load_dataset("splat.sim")))
    # Every gene but the two without counts was fitted
    fitted <- setdiff(genes, c("Gene1", "Gene2"))
    expect_identical(names(ob@Profile@p_values), fitted)
    expect_identical(res[fitted, "p_value"], unname(ob@Profile@p_values))
    expect_identical(
        res[fitted, "adj_p_value"], unname(ob@Profile@adj_p_values)
    )
    expect_identical(res$non_flat, genes %in% ob@Profile@non_flat)
    sol <- showSol(ob)
    expect_identical(res[rownames(sol), "model_p_value"], sol[["p-value"]])
    expect_identical(res[rownames(sol), "r_squared"], sol[["R-squared"]])
    expect_true(all(is.na(res[!genes %in% rownames(sol), "r_squared"])))
    expect_true(any(res$influential))
    expect_identical(res$influential, genes %in% colnames(showInflu(ob)))
    for (name in names(ob@Significant@genes)) {
        expect_identical(
            res[[paste0("sig_", name)]], genes %in% ob@Significant@genes[[name]]
        )
    }
    clusters <- ob@Significant@clusters
    expect_type(res$cluster, "integer")
    expect_identical(
        res[names(clusters), "cluster"], unname(unlist(clusters))
    )
    expect_true(all(is.na(res[!genes %in% names(clusters), "cluster"])))
    # The genes left out by sc.p.vector(): NA values, FALSE memberships
    zero <- res[c("Gene1", "Gene2"), ]
    for (col in c("p_value", "adj_p_value", "model_p_value", "r_squared")) {
        expect_identical(zero[[col]], c(NA_real_, NA_real_))
    }
    expect_identical(zero$cluster, c(NA_integer_, NA_integer_))
    logical_cols <- colnames(res)[vapply(as.list(res), is.logical, NA)]
    expect_length(logical_cols, 2 + length(ob@Significant@genes))
    for (col in logical_cols) {
        expect_false(anyNA(res[[col]]))
        expect_identical(zero[[col]], c(FALSE, FALSE))
    }
})

test_that("scmpToRowData() writes scmpResults() into rowData() by name", {
    ob <- splat_clustered()
    res <- scmpResults(ob)
    genes <- rownames(res)
    # The genes in reverse order, five genes of the object missing and two
    # genes that the object does not have; a column of the user and a
    # column left by an earlier call
    kept <- rev(genes[-(1:5)])
    rows <- c(kept[1:10], "NewGene1", kept[-(1:10)], "NewGene2")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1,
            nrow = length(rows), ncol = 2,
            dimnames = list(rows, c("a", "b"))
        )),
        rowData = S4Vectors::DataFrame(symbol = tolower(rows), scmp_old = 1)
    )
    out <- scmpToRowData(sce, ob)
    rd <- SummarizedExperiment::rowData(out)
    expect_identical(rownames(out), rows)
    expect_identical(
        colnames(rd), c("symbol", paste0("scmp_", colnames(res)))
    )
    expect_identical(rd$symbol, tolower(rows))
    # Read back: the genes of the object equal scmpResults(), the others
    # are NA in every column
    back <- rd[, -1]
    colnames(back) <- sub("^scmp_", "", colnames(back))
    in_ob <- rows %in% genes
    expect_identical(
        as.data.frame(back[in_ob, ]), as.data.frame(res[rows[in_ob], ])
    )
    expect_true(all(is.na(as.data.frame(back[!in_ob, ]))))
    expect_identical(
        SummarizedExperiment::assay(out), SummarizedExperiment::assay(sce)
    )
    # A later call replaces the columns: without results, none is left
    out2 <- scmpToRowData(out, splat_scmp())
    expect_identical(
        colnames(SummarizedExperiment::rowData(out2)), "symbol"
    )
})

test_that("scmpToRowData() works on the SCE the analysis started from", {
    sce <- load_dataset("splat.sim")
    ob <- splat_filtered()
    out <- scmpToRowData(sce, ob)
    rd <- SummarizedExperiment::rowData(out)
    old <- colnames(SummarizedExperiment::rowData(sce))
    expect_identical(
        colnames(rd), c(old, paste0("scmp_", colnames(scmpResults(ob))))
    )
    expect_identical(rd$scmp_p_value, scmpResults(ob)$p_value)
    expect_identical(rd$Gene, SummarizedExperiment::rowData(sce)$Gene)
})

test_that("scmpResults() and scmpToRowData() check their arguments", {
    msg_class <- "Please provide object of class 'ScMaSigPro'"
    msg_sce <- paste(
        "'sce' must be a SingleCellExperiment or another",
        "SummarizedExperiment."
    )
    msg_none <- paste(
        "None of the genes of 'sce' (its row names) is a gene of",
        "'scmpObj'."
    )
    sce <- load_dataset("splat.sim")
    expect_error(scmpResults(sce), msg_class, fixed = TRUE)
    expect_error(scmpToRowData(sce, sce), msg_class, fixed = TRUE)
    expect_error(
        scmpToRowData(SummarizedExperiment::assay(sce), splat_filtered()),
        msg_sce,
        fixed = TRUE
    )
    other <- sce
    rownames(other) <- paste0("X", seq_len(nrow(sce)))
    expect_error(scmpToRowData(other, splat_filtered()), msg_none,
        fixed = TRUE
    )
    rownames(other) <- NULL
    expect_error(scmpToRowData(other, splat_filtered()), msg_none,
        fixed = TRUE
    )
})
