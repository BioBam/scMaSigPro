# sc.squeeze() must reproduce a manual pseudo-bulk: for every bin, the sum (or
# rounded mean) of the counts of its member cells, for every gene.

bin_members <- function(ob) {
    bins <- cDense(ob)
    members <- strsplit(bins[[ob@Parameters@bin_mem_col]], "|", fixed = TRUE)
    names(members) <- rownames(bins)
    members
}

test_that("pseudo-bulk counts equal the per-bin sums of the member cells", {
    ob <- splat_binned()
    bins <- cDense(ob)
    counts <- eSparse(ob)
    bulk <- eDense(ob)
    members <- bin_members(ob)

    # every cell is in exactly one bin, and bin sizes are recorded correctly
    expect_setequal(unlist(members, use.names = FALSE), colnames(counts))
    expect_identical(anyDuplicated(unlist(members, use.names = FALSE)), 0L)
    expect_equal(unname(lengths(members)), as.numeric(bins[[ob@Parameters@bin_size_col]]))

    manual <- vapply(members, function(cells) {
        rowSums(counts[, cells, drop = FALSE])
    }, numeric(nrow(counts)))
    expect_equal(manual, bulk[, names(members)])
    expect_equal(sum(bulk), sum(counts))
})

test_that("aggregate = 'mean' gives the rounded per-bin mean of the member cells", {
    ob <- sc.squeeze(splat_scmp(), aggregate = "mean")
    counts <- eSparse(ob)
    bulk <- eDense(ob)
    members <- bin_members(ob)

    manual <- vapply(members, function(cells) {
        round(rowMeans(counts[, cells, drop = FALSE]))
    }, numeric(nrow(counts)))
    expect_equal(manual, bulk[, names(members)])
})

test_that("pb_counts() can be re-run on an object and gives the same matrix", {
    ob <- splat_binned()
    again <- pb_counts(ob, assay_name = "counts", cluster_count_by = "sum")
    expect_equal(eDense(again), eDense(ob))
    expect_equal(cDense(again), cDense(ob))
})
