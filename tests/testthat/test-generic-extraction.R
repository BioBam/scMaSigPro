# Accessor generics: getters return the expected matrices and data frames,
# setters round-trip.

test_that("eSparse()/eDense() return the count matrices with dimnames", {
    ob <- splat_binned()
    raw <- as.matrix(SingleCellExperiment::counts(load_dataset("splat.sim")))

    expect_true(is.matrix(eSparse(ob)))
    expect_equal(eSparse(ob)[, colnames(raw)], raw)

    bulk <- eDense(ob)
    expect_true(is.matrix(bulk))
    expect_identical(dim(bulk), c(nrow(raw), nrow(cDense(ob))))
    expect_identical(rownames(bulk), rownames(raw))
    expect_identical(colnames(bulk), rownames(cDense(ob)))
})

test_that("named assays can be set and retrieved in both slots", {
    ob <- splat_binned()

    m <- log1p(eSparse(ob))
    expect_message(eSparse(ob, "logcounts") <- m, "Creating a new assay")
    expect_equal(eSparse(ob, "logcounts"), m)
    expect_warning(eSparse(ob, "logcounts") <- m, "Overwriting")

    d <- eDense(ob) * 2
    expect_message(eDense(ob, "scaled") <- d, "Creating a new assay")
    expect_equal(eDense(ob, "scaled"), d)
    expect_equal(eDense(ob), eDense(splat_binned()))
})

test_that("cSparse()/cDense() return data frames keyed by cell and bin", {
    ob <- splat_binned()

    cells <- cSparse(ob)
    expect_s3_class(cells, "data.frame")
    expect_true(all(c("Pseudotime", "Path", "scmp_binned_pseudotime") %in% names(cells)))
    expect_setequal(rownames(cells), colnames(eSparse(ob)))

    bins <- cDense(ob)
    expect_s3_class(bins, "data.frame")
    expect_true(all(c(
        ob@Parameters@bin_col, ob@Parameters@bin_mem_col, ob@Parameters@bin_size_col,
        ob@Parameters@bin_ptime_col, ob@Parameters@path_col
    ) %in% names(bins)))
    expect_identical(rownames(bins), bins[[ob@Parameters@bin_col]])
})

test_that("cDense<- and cSparse<- replace the metadata", {
    ob <- splat_binned()

    bins <- cDense(ob)
    bins$extra <- seq_len(nrow(bins))
    cDense(ob) <- bins
    expect_identical(cDense(ob)$extra, seq_len(nrow(bins)))

    cells <- cSparse(ob)
    cells$flag <- TRUE
    cSparse(ob) <- cells
    expect_true(all(cSparse(ob)$flag))
})

test_that("predictors()/pathAssign() expose the design after sc.set.poly()", {
    ob <- splat_design(2)
    expect_true(is.matrix(predictors(ob)))
    expect_identical(ncol(predictors(ob)), 5L)
    expect_identical(rownames(predictors(ob)), rownames(cDense(ob)))
    expect_true(is.matrix(pathAssign(ob)))
    expect_identical(nrow(pathAssign(ob)), nrow(cDense(ob)))

    p <- predictors(ob)
    predictors(ob) <- p * 1
    expect_equal(predictors(ob), p)
})

test_that("cSparse<- rejects a cell table that does not match the cells", {
    ob <- splat_binned()
    cells <- cSparse(ob)
    expect_error(
        cSparse(ob) <- cells[-1, ],
        paste(
            "'value' must have one row per cell of the 'Sparse' slot:",
            "199 rows for 200 cells."
        ),
        fixed = TRUE
    )
    expect_error(cSparse(ob, cells[-1, ]), "one row per cell", fixed = TRUE)
    msg <- paste(
        "The rownames of 'value' must be the cell names of the 'Sparse'",
        "slot, in the same order."
    )
    expect_error(cSparse(ob) <- cells[rev(rownames(cells)), ], msg,
        fixed = TRUE
    )
    rownames(cells) <- paste0("x", rownames(cells))
    expect_error(cSparse(ob) <- cells, msg, fixed = TRUE)
    expect_identical(ob, splat_binned())
})

test_that("cSparse<- takes a table without rownames in cell order", {
    ob <- splat_binned()
    cells <- cSparse(ob)
    unnamed <- cells
    rownames(unnamed) <- NULL
    cSparse(ob) <- unnamed
    expect_identical(cSparse(ob), cells)
    expect_sparse_aligned(ob)
})

test_that("eSparse<- and eDense<- reject a matrix that does not fit the slot", {
    ob <- splat_binned()
    m <- eSparse(ob)
    expect_error(
        eSparse(ob, "x") <- m[, -1],
        paste(
            "'value' must have one row per gene and one column per cell of",
            "the 'Sparse' slot: 100 x 199 given, 100 x 200 needed."
        ),
        fixed = TRUE
    )
    expect_error(eSparse(ob, "x") <- t(m), "one column per cell", fixed = TRUE)
    expect_error(
        eSparse(ob, "x") <- m[rev(rownames(m)), ],
        paste(
            "The rownames of 'value' must be the gene names of the 'Sparse'",
            "slot, in the same order."
        ),
        fixed = TRUE
    )
    expect_error(
        eSparse(ob, "counts") <- m[, rev(colnames(m))],
        paste(
            "The colnames of 'value' must be the cell names of the 'Sparse'",
            "slot, in the same order."
        ),
        fixed = TRUE
    )
    d <- eDense(ob)
    expect_error(
        eDense(ob, "x") <- d[-1, ],
        paste(
            "'value' must have one row per gene and one column per bin of",
            "the 'Dense' slot: 99 x 16 given, 100 x 16 needed."
        ),
        fixed = TRUE
    )
    expect_error(
        eDense(ob, "bulk.counts") <- d[, rev(colnames(d))],
        paste(
            "The colnames of 'value' must be the bin names of the 'Dense'",
            "slot, in the same order."
        ),
        fixed = TRUE
    )
    expect_identical(ob, splat_binned())
})

test_that("eSparse<- and eDense<- fill in missing dimnames", {
    ob <- splat_binned()
    m <- log1p(eSparse(ob))
    expect_message(eSparse(ob, "logcounts") <- unname(m), "Creating")
    expect_identical(eSparse(ob, "logcounts"), m)
    d <- eDense(ob) * 2
    expect_message(eDense(ob, "scaled") <- unname(d), "Creating")
    expect_identical(eDense(ob, "scaled"), d)
})

test_that("cDense<- with another number of bins drops the counts", {
    # As when bins are rebuilt by hand (the Shared root vignette): the
    # object stays valid and pb_counts() computes the counts of the new bins.
    ob <- splat_binned()
    bins <- cDense(ob)[-1, ]
    expect_message(
        cDense(ob) <- bins,
        paste(
            "'value' has 15 bins instead of 16: the pseudo-bulk counts of the",
            "'Dense' slot are removed; run pb_counts() to compute them for",
            "the new bins."
        ),
        fixed = TRUE
    )
    expect_valid_scmp(ob)
    expect_identical(cDense(ob), bins)
    expect_error(eDense(ob), "No assays found in the object.", fixed = TRUE)
    ob <- pb_counts(ob)
    expect_identical(eDense(ob), eDense(splat_binned())[, rownames(bins)])
})

test_that("pathAssign<- and predictors<- reject a mismatched matrix", {
    ob <- splat_design(2)
    p <- predictors(ob)
    pa <- pathAssign(ob)
    rows <- paste(
        "slot 'Design': 'assignment_matrix' must have the rows of",
        "'predictor_matrix'"
    )
    expect_error(predictors(ob) <- p[-1, ], rows, fixed = TRUE)
    expect_error(predictors(ob) <- unname(p), rows, fixed = TRUE)
    expect_error(predictors(ob, p[-1, ]), rows, fixed = TRUE)
    expect_error(pathAssign(ob) <- pa[-1, ], rows, fixed = TRUE)
    expect_error(pathAssign(ob, pa[-1, ]), rows, fixed = TRUE)
    expect_error(
        predictors(ob) <- p[, -1],
        paste(
            "slot 'Design': 'groups.vector' must have one value per column",
            "of 'predictor_matrix'"
        ),
        fixed = TRUE
    )
    expect_identical(ob, splat_design(2))
})

test_that("every setter refuses to return an invalid object", {
    # The setters check the whole object, so a problem in another slot is
    # reported too.
    ob <- splat_design(2)
    ob@Parameters@anno_col <- character()
    msg <- "slot 'Parameters': 'anno_col' must be a single non-missing value"
    expect_error(cSparse(ob) <- cSparse(ob), msg, fixed = TRUE)
    expect_error(cSparse(ob, cSparse(ob)), msg, fixed = TRUE)
    expect_error(cDense(ob) <- cDense(ob), msg, fixed = TRUE)
    expect_error(cDense(ob, cDense(ob)), msg, fixed = TRUE)
    expect_error(
        suppressMessages(eSparse(ob, "x") <- eSparse(ob)), msg,
        fixed = TRUE
    )
    expect_error(
        suppressMessages(eDense(ob, "x") <- eDense(ob)), msg,
        fixed = TRUE
    )
    expect_error(pathAssign(ob) <- pathAssign(ob), msg, fixed = TRUE)
    expect_error(pathAssign(ob, pathAssign(ob)), msg, fixed = TRUE)
    expect_error(predictors(ob) <- predictors(ob), msg, fixed = TRUE)
    expect_error(predictors(ob, predictors(ob)), msg, fixed = TRUE)
})
