# create_scmp() with bins made outside scMaSigPro (roadmap #18d, D18)

# The cells of splat.sim and their bins, as a user would bring them
own_bins <- function(scmp, binned) {
    bins <- cDense(binned)
    list(
        counts = eSparse(scmp), cell_data = cSparse(scmp),
        bin_counts = eDense(binned),
        bin_cell_data = data.frame(
            Pseudotime = bins$scmp_binned_pseudotime, Path = bins$Path,
            scmp_bin_size = bins$scmp_bin_size, row.names = rownames(bins)
        )
    )
}

test_that("create_scmp() builds the Dense slot from given bins (ledger q10)", {
    x <- own_bins(splat_scmp(), splat_binned())
    ob <- create_scmp(
        counts = x$counts, cell_data = x$cell_data,
        bin_counts = x$bin_counts, bin_cell_data = x$bin_cell_data,
        ptime_col = "Pseudotime", path_col = "Path"
    )
    expect_true(validObject(ob))
    expect_identical(eDense(ob), x$bin_counts)
    expect_identical(ob@Parameters@bin_ptime_col, "Pseudotime")
    # The workflow gives the results of the bins of sc.squeeze()
    ob <- sc.p.vector(sc.set.poly(ob, poly_degree = 2), verbose = FALSE)
    ob <- sc.t.fit(ob, verbose = FALSE)
    ref <- splat_fitted("default")
    expect_identical(ob@Profile, ref@Profile)
    for (s in c(
        "significance_matrix", "coefficient_matrix", "t_score_matrix",
        "path_coefficient_matrix", "influential"
    )) {
        expect_identical(unname(slot(ob@Estimate, s)),
            unname(slot(ref@Estimate, s)),
            info = s
        )
    }
    # Without bin sizes, the sizes are unknown
    x$bin_cell_data$scmp_bin_size <- NULL
    ob <- create_scmp(
        counts = x$counts, cell_data = x$cell_data,
        bin_counts = x$bin_counts, bin_cell_data = x$bin_cell_data,
        ptime_col = "Pseudotime", path_col = "Path"
    )
    expect_true(all(is.na(cDense(ob)$scmp_bin_size)))
})

test_that("create_scmp() checks the bins it is given", {
    x <- own_bins(splat_scmp(), splat_binned())
    make <- function(...) {
        args <- list(
            counts = x$counts, cell_data = x$cell_data,
            bin_counts = x$bin_counts, bin_cell_data = x$bin_cell_data,
            ptime_col = "Pseudotime", path_col = "Path"
        )
        new <- list(...)
        args[names(new)] <- new
        do.call(create_scmp, args)
    }
    expect_error(make(bin_cell_data = NULL),
        "Please provide both 'bin_counts' and 'bin_cell_data'.",
        fixed = TRUE
    )
    expect_error(make(use_as_bin = TRUE),
        "'bin_counts' cannot be combined with 'use_as_bin = TRUE'.",
        fixed = TRUE
    )
    expect_error(make(bin_counts = x$bin_counts[, -1]),
        "Number of bins in bin_counts and bin_cell_data are different.",
        fixed = TRUE
    )
    expect_error(make(bin_counts = x$bin_counts[, rev(colnames(x$bin_counts))]),
        "Colnames of bin_counts and rownames of bin_cell_data are different.",
        fixed = TRUE
    )
    expect_error(make(bin_counts = x$bin_counts[rev(rownames(x$bin_counts)), ]),
        "The rows of bin_counts must be the genes of counts.",
        fixed = TRUE
    )
    expect_error(make(bin_cell_data = x$bin_cell_data[, -2]),
        "Column name 'Path' not found in bin_cell_data.",
        fixed = TRUE
    )
})

test_that("create_scmp() stores feature_data as rowData (ledger q10)", {
    x <- own_bins(splat_scmp(), splat_binned())
    genes <- data.frame(symbol = tolower(rownames(x$counts)))
    ob <- create_scmp(
        counts = x$counts, cell_data = x$cell_data, feature_data = genes,
        bin_counts = x$bin_counts, bin_cell_data = x$bin_cell_data,
        ptime_col = "Pseudotime", path_col = "Path"
    )
    expect_identical(
        SummarizedExperiment::rowData(ob@Sparse)$symbol, genes$symbol
    )
    expect_identical(
        rownames(SummarizedExperiment::rowData(ob@Dense)), rownames(x$counts)
    )
    expect_error(
        create_scmp(
            counts = x$counts, cell_data = x$cell_data,
            feature_data = genes[-1, , drop = FALSE],
            ptime_col = "Pseudotime", path_col = "Path"
        ),
        "feature_data must have one row per gene of counts.",
        fixed = TRUE
    )
})
