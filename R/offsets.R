#' Size factors by the median-of-ratios method
#'
#' Estimates one size factor per sample (column) of a count matrix with the
#' median-of-ratios method of Anders and Huber (2010), written from the
#' published definition. For the count `k_ij` of gene `i` in sample `j`, the
#' size factor of sample `j` is the median over genes of `k_ij / g_i`, where
#' `g_i` is the geometric mean of gene `i` over the `m` samples. The
#' computation is on the log scale, `log(g_i) = mean_j(log(k_ij))` and
#' `s_j = exp(median_i(log(k_ij) - log(g_i)))`, so with an even number of
#' genes the median is the geometric mean of the two middle ratios.
#'
#' * `type = "ratio"`: `g_i = (k_i1 * ... * k_im)^(1/m)`. A gene with a zero
#'   or missing count in any sample has no finite positive geometric mean
#'   and is skipped.
#' * `type = "poscounts"` (the positive-counts variant of the DESeq2 package,
#'   Love et al. 2014): `g_i` is the `m`-th root of the product of the
#'   positive counts of gene `i`, so genes with some zeros are used, and the
#'   median for sample `j` is taken over the genes with a positive count in
#'   that sample. Genes whose counts are all zero, or that have a missing
#'   count, are skipped. A sample with no positive count among the usable
#'   genes gets `NA`.
#'
#' The input is converted with [as.matrix()] first, so a base matrix, a data
#' frame or a `Matrix` object (such as the `dgeMatrix` that `dgCMatrix + 1`
#' returns) gives the same result.
#'
#' @param counts A matrix-like object of non-negative counts, genes in rows
#'   and samples in columns.
#' @param type `"ratio"` (default) or `"poscounts"`, as described above.
#'
#' @return A numeric vector with one size factor per column of `counts`,
#'   named by its column names. It stops with an error if `counts` is not
#'   two-dimensional or if no gene can be used.
#'
#' @references Anders, S. and Huber, W. (2010). Differential expression
#' analysis for sequence count data. Genome Biology 11, R106.
#' \doi{10.1186/gb-2010-11-10-r106}
#'
#' Love, M.I., Huber, W. and Anders, S. (2014). Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550.
#' \doi{10.1186/s13059-014-0550-8}
#'
#' @importFrom stats median
#' @noRd
.size_factors <- function(counts, type = c("ratio", "poscounts")) {
    type <- match.arg(type)
    # A plain vector stops with the error base rowMeans() gave
    .check(length(dim(counts)) == 2L,
        msg = "'x' must be an array of at least two dimensions"
    )
    counts <- as.matrix(counts)
    log_counts <- log(counts)

    # Log geometric mean of every gene over all samples
    if (type == "ratio") {
        # Why: log(0) is -Inf, so a gene with a zero count has no finite
        # geometric mean and stays out of the medians.
        log_geo_means <- rowMeans(log_counts)
    } else {
        log_pos <- log_counts
        log_pos[counts <= 0] <- 0
        log_geo_means <- rowMeans(log_pos)
        log_geo_means[rowSums(counts > 0) == 0] <- -Inf
    }
    use <- is.finite(log_geo_means)
    .check(any(use), msg = if (type == "ratio") {
        paste(
            "Cannot compute size factors: every gene contains at least one",
            "zero or missing count."
        )
    } else {
        paste(
            "Cannot compute size factors: every gene has only zero counts",
            "or a missing count."
        )
    })

    # Median over the usable genes of the log ratio to the geometric mean
    sf <- vapply(seq_len(ncol(counts)), function(j) {
        keep <- use & counts[, j] > 0
        exp(median(log_counts[keep, j] - log_geo_means[keep]))
    }, numeric(1))
    names(sf) <- colnames(counts)
    sf
}

# The offset of sc.p.vector() and sc.t.fit() (Srivastava et al. (2024),
# Section 2.3 and Supplementary Section S1.3): the logarithm of the size
# factor of each bin, from the raw counts of the genes that sc.p.vector()
# fits. Only the genes with a count above 0 in every bin have a geometric
# mean above 0, so at least one is needed.
.bin_offset <- function(counts) {
    counts <- as.matrix(counts)
    .check(isTRUE(any(rowSums(!(counts > 0)) == 0)),
        msg = paste(
            "Cannot compute the size factors of the offset: every gene has a",
            "zero or missing count in at least one bin of the design. Fewer",
            "and larger bins (for example sc.squeeze() with a smaller",
            "'drop_fac') can leave genes without a zero count; otherwise fit",
            "without the offset, with offset = FALSE. 'min_counts' and",
            "'min_bins_expressed' do not help: the size factors are computed",
            "before them."
        )
    )
    # Why: the offset is added to the log of the mean.
    log(.size_factors(counts))
}
