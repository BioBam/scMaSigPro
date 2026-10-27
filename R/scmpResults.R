#' @title Per-gene results of an analysis
#'
#' @description
#' `scmpResults()` collects the per-gene results stored in a
#' [ScMaSigPro-class] object in one table with a row per gene.
#' `scmpToRowData()` writes that table into the `rowData()` of a
#' `SingleCellExperiment`, such as the one the analysis started from, so that
#' the results travel with the data.
#'
#' @param sce A `SingleCellExperiment`, or another `SummarizedExperiment`,
#'   whose row names are gene names.
#'
#' @details
#' `scmpResults()` has one row per gene of the object (the rows of the
#' `Sparse` slot), named by gene, and the columns of the steps that have been
#' run:
#' * [sc.p.vector()]: `p_value` and `adj_p_value`, the p-value of the full
#'   model and its adjusted value, `NA` for the genes that were not fitted
#'   (too few values or only zero counts); `non_flat`, `TRUE` for the genes
#'   it selected.
#' * [sc.t.fit()]: `model_p_value` and `r_squared`, the p-value and the
#'   R-squared of the model selected by the stepwise regression (the
#'   `p-value` and `R-squared` columns of [showSol()]), `NA` for the genes
#'   without a row there; `influential`, `TRUE` for the genes with
#'   influential observations (the columns of [showInflu()]).
#' * [sc.filter()]: one column per gene list, named `sig_` followed by the
#'   name of the list, `TRUE` for the genes of the list.
#' * [sc.cluster.trend()]: `cluster`, the cluster of each clustered gene, `NA`
#'   for the other genes.
#'
#' The logical columns are `FALSE`, never `NA`, for the genes outside a set.
#'
#' `scmpToRowData()` matches the genes by name. It removes the columns of
#' `rowData(sce)` whose names start with `scmp_`, as an earlier call writes
#' them, and adds the columns of `scmpResults()` with that prefix, such as
#' `scmp_p_value`. The genes of `sce` that are not in the object get `NA` in
#' every new column, and the genes of the object that are not in `sce` are
#' left out. It stops if no gene matches.
#'
#' @return `scmpResults()`: a
#' \code{\link[S4Vectors:DataFrame-class]{DataFrame}}. `scmpToRowData()`:
#' `sce` with the new `rowData()` columns.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' res <- scmpResults(scmp.ob)
#' head(res)
#' # The genes of the "Path1" list, ordered by R-squared
#' sig <- res[res$sig_Path1, ]
#' head(sig[order(sig$r_squared, decreasing = TRUE), ])
#'
#' # scmp.ob was made from the bundled splat.sim data set
#' data("splat.sim", package = "scMaSigPro")
#' sce <- scmpToRowData(splat.sim, scmp.ob)
#' head(SummarizedExperiment::rowData(sce)[, c("Gene", "scmp_p_value")])
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @family results
#' @inheritParams scmp-params
#'
#' @export
scmpResults <- function(scmpObj) {
    .check_scmp(scmpObj)
    genes <- rownames(scmpObj@Sparse)
    .check(!is.null(genes) || nrow(scmpObj@Sparse) == 0L,
        msg = "The genes of the 'Sparse' slot have no names."
    )
    res <- DataFrame(row.names = genes)

    # sc.p.vector(): the full model of every fitted gene
    profile <- scmpObj@Profile
    if (length(profile@p_values) > 0L) {
        res[["p_value"]] <- unname(profile@p_values[genes])
        res[["adj_p_value"]] <- unname(profile@adj_p_values[genes])
        res[["non_flat"]] <- genes %in% profile@non_flat
    }

    # sc.t.fit(): the selected model of every gene in sol
    sol <- scmpObj@Estimate@significance_matrix
    if (nrow(sol) > 0L) {
        in_sol <- match(genes, rownames(sol))
        res[["model_p_value"]] <- unname(sol[in_sol, "p-value"])
        res[["r_squared"]] <- unname(sol[in_sol, "R-squared"])
        res[["influential"]] <- genes %in%
            colnames(scmpObj@Estimate@influential)
    }

    # sc.filter() and sc.cluster.trend(): gene lists and clusters
    gene_lists <- scmpObj@Significant@genes
    for (name in names(gene_lists)) {
        res[[paste0("sig_", name)]] <- genes %in% gene_lists[[name]]
    }
    clusters <- scmpObj@Significant@clusters
    if (length(clusters) > 0L) {
        cluster <- as.integer(unlist(clusters, use.names = FALSE))
        res[["cluster"]] <- cluster[match(genes, names(clusters))]
    }
    res
}

#' @rdname scmpResults
#' @export
scmpToRowData <- function(sce, scmpObj) {
    .check_scmp(scmpObj)
    .check(is(sce, "SummarizedExperiment"),
        msg = paste(
            "'sce' must be a SingleCellExperiment or another",
            "SummarizedExperiment."
        )
    )
    res <- scmpResults(scmpObj)
    genes <- rownames(sce)
    in_res <- match(genes, rownames(res))
    .check(!is.null(genes) && any(!is.na(in_res)),
        msg = paste(
            "None of the genes of 'sce' (its row names) is a gene of",
            "'scmpObj'."
        )
    )

    # Replace the columns of an earlier call
    row_data <- rowData(sce)
    row_data <- row_data[, !startsWith(colnames(row_data), "scmp_"),
        drop = FALSE
    ]
    for (name in colnames(res)) {
        row_data[[paste0("scmp_", name)]] <- res[[name]][in_res]
    }
    rowData(sce) <- row_data
    sce
}
