# Structural validity of the object and of its two SingleCellExperiment slots.
# Regression tests for roadmap items S1 and S2 (invalid Dense slot, cell
# metadata misaligned with the count matrix after binning).

test_that("objects are valid S4 at every stage of the splat.sim workflow", {
    expect_valid_scmp(splat_scmp())
    expect_valid_scmp(splat_binned())
    expect_valid_scmp(splat_design(2))
    expect_valid_scmp(splat_fitted("default"))
    expect_valid_scmp(splat_filtered())
    expect_valid_scmp(splat_clustered())
})

test_that("objects are valid S4 after sc.squeeze() and sc.restruct() on multi.lin.sce", {
    expect_valid_scmp(multi_binned())
    expect_valid_scmp(multi_restruct())
})

test_that("the bundled scmp.ob is a valid object", {
    expect_valid_scmp(load_dataset("scmp.ob"))
})

test_that("cell metadata stays aligned with the counts after binning", {
    ob <- splat_binned()
    expect_sparse_aligned(ob)

    raw <- as.matrix(SingleCellExperiment::counts(load_dataset("splat.sim")))
    via_accessor <- as.matrix(SingleCellExperiment::counts(ob@Sparse))
    expect_equal(via_accessor[, colnames(raw)], raw)
    expect_equal(eSparse(ob)[, colnames(raw)], raw)

    meta <- cSparse(ob)
    expect_identical(rownames(meta), colnames(via_accessor))
})

test_that("the Dense slot has gene rows and bin columns", {
    for (ob in list(splat_binned(), multi_binned(), multi_restruct())) {
        expect_identical(dim(ob@Dense), dim(eDense(ob)))
        expect_identical(rownames(ob@Dense), rownames(eDense(ob)))
        expect_identical(colnames(ob@Dense), rownames(cDense(ob)))
        expect_identical(colnames(eDense(ob)), rownames(cDense(ob)))
    }
})

test_that("sc.squeeze() realigns cell metadata that arrives out of order", {
    # Objects binned by 0.0.4 carry cell metadata in a different order from the
    # counts; binning such an object again must repair the pairing.
    ob <- splat_scmp()
    ob@Sparse@colData <- ob@Sparse@colData[rev(seq_len(ncol(ob@Sparse))), ]
    ob <- sc.squeeze(ob)
    expect_valid_scmp(ob)
    expect_sparse_aligned(ob)
    expect_equal(eDense(ob), eDense(splat_binned()))
})

test_that("pb_counts() names the bins after the rownames of the Dense metadata", {
    ob <- splat_binned()
    bins <- cDense(ob)
    rownames(bins) <- paste0("renamed_", rownames(bins))
    cDense(ob) <- bins
    ob <- pb_counts(ob)
    expect_valid_scmp(ob)
    expect_identical(colnames(eDense(ob)), rownames(cDense(ob)))
    expect_equal(unname(eDense(ob)), unname(eDense(splat_binned())))
    ob <- sc.set.poly(ob, poly_degree = 2)
    expect_identical(rownames(predictors(ob)), colnames(eDense(ob)))
    expect_no_error(sc.p.vector(ob, verbose = FALSE))
})

test_that("pb_counts() falls back to the bin column when bins have no rownames", {
    ob <- splat_binned()
    bins <- cDense(ob)
    rownames(bins) <- NULL
    cDense(ob) <- bins
    ob <- pb_counts(ob)
    expect_valid_scmp(ob)
    expect_identical(colnames(eDense(ob)), cDense(ob)[[ob@Parameters@bin_col]])
})

test_that("pb_counts() rejects duplicated bin names", {
    ob <- splat_binned()
    bins <- SummarizedExperiment::colData(ob@Dense)
    rownames(bins)[2] <- rownames(bins)[1]
    ob@Dense@colData <- bins
    expect_error(pb_counts(ob), "must be unique")
})

test_that("sc.squeeze() removes cells that fall outside every bin, with a warning", {
    # split_bins leaves gaps between the split bins on multi.lin.sce (3 cells
    # with Sturges); 0.0.4 dropped them from the metadata only, leaving an
    # invalid object.
    expect_warning(
        ob <- sc.squeeze(multi_scmp(), split_bins = TRUE),
        "could not be assigned to a bin"
    )
    expect_valid_scmp(ob)
    expect_sparse_aligned(ob)
    expect_identical(ncol(ob@Sparse), 1974L)
})
