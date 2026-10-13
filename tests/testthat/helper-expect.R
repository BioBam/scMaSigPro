# Custom expectations shared by the test files.

#' Expect that an ScMaSigPro object and its two SingleCellExperiment slots are
#' valid S4 objects.
expect_valid_scmp <- function(object) {
    act <- testthat::quasi_label(rlang::enquo(object), arg = "object")
    msg <- tryCatch(
        {
            methods::validObject(act$val)
            methods::validObject(act$val@Sparse)
            methods::validObject(act$val@Dense)
            NULL
        },
        error = function(e) conditionMessage(e)
    )
    testthat::expect(
        is.null(msg),
        sprintf("%s is not a valid ScMaSigPro object: %s", act$lab, msg)
    )
    invisible(act$val)
}

#' Expect that the cell-level metadata of the Sparse slot is aligned, cell by
#' cell, with the columns of the count matrix.
expect_sparse_aligned <- function(object) {
    act <- testthat::quasi_label(rlang::enquo(object), arg = "object")
    meta_cells <- rownames(SummarizedExperiment::colData(act$val@Sparse))
    assay_cells <- colnames(SummarizedExperiment::assay(act$val@Sparse, withDimnames = FALSE))
    testthat::expect(
        identical(meta_cells, assay_cells),
        sprintf(
            "%s: colData rows are not in the same order as the assay columns (first mismatch at %d)",
            act$lab, which(meta_cells != assay_cells)[1]
        )
    )
    invisible(act$val)
}
