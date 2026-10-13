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
