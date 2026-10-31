#' @title Update an ScMaSigPro object to the current class definition
#'
#' @description
#' `updateObject()` brings a [ScMaSigPro-class] object saved by an earlier
#' version of scMaSigPro, for example with `saveRDS()` or `save()`, up to
#' date. Call it on such an object after `readRDS()` or `load()`, before the
#' next step of the workflow.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param ... Passed to the `updateObject()` methods of the `Sparse` and
#'   `Dense` slots.
#' @param verbose Whether to print a message when the object is updated.
#'
#' @details
#' Objects saved before the `version` slot was added (scMaSigPro 0.0.5 and
#' earlier) fail `validObject()`; `updateObject()` adds the slot with the
#' version of the package that is running. For every object, the `Sparse` and
#' `Dense` slots are passed to the `updateObject()` method of
#' SingleCellExperiment, and the result is checked with `validObject()`. A
#' current object, such as one created by this version of the package, is
#' returned unchanged.
#'
#' Objects saved before the `min_counts` and `min_bins_expressed` parameters
#' existed get both set to 0, the values that keep every gene
#' ([sc.p.vector()] had no such filter then), and the version of the package
#' that is running.
#'
#' Objects saved before the `test` parameter existed get it set to `"chisq"`,
#' the test of [sc.p.vector()] then, and the version of the package that is
#' running.
#'
#' Objects saved with the `log_offset` parameter lose it. The fits of an
#' object saved with `log_offset = FALSE`, the default, used the size factor
#' of each bin as the offset: `updateObject()` stores the logarithm of the
#' size factors instead, the offset of [sc.p.vector()] and [sc.t.fit()], and
#' warns that the fits in the `Profile` and `Estimate` slots were made with
#' the size factors, until [sc.p.vector()] is run again. The offsets of an
#' object saved with `log_offset = TRUE` are kept. An object whose fits used
#' no offset, or that was not fitted, keeps no offset: the offsets that an
#' earlier fit left in it, or the zeros of [sc.set.poly()], are removed.
#'
#' Objects fitted with `link = "identity"` by scMaSigPro 0.0.5 or earlier
#' hold a family labelled with the identity link but with the link functions
#' of the log link, with which their fits were made: `updateObject()` warns,
#' and [sc.t.fit()] stops on them, until [sc.p.vector()] is run again.
#'
#' Objects binned by scMaSigPro 0.0.4 cannot be updated: their `Dense` slot
#' is not a valid `SingleCellExperiment`, and the cell metadata of their
#' `Sparse` slot may not be in the order of the counts. `updateObject()`
#' stops with the command that bins the cells of the `Sparse` slot again with
#' [sc.squeeze()]; the later steps of the workflow must then be run again.
#'
#' Earlier versions kept the results of the later steps when a step was run
#' again, so an object saved after that may hold results computed from the
#' slots that the step replaced. `updateObject()` warns, without changing
#' the object, and names the step to run again, when the design does not
#' have the binned pseudotime of the bins of the `Dense` slot, when tests
#' made with an offset find none in the design, when the models are not
#' those of non-flat genes with a coefficient for each column of the
#' design, when the gene lists hold genes without a model, or when the
#' clusters hold genes that are not in the lists. It cannot tell every such
#' result from a current one, such as tests of counts that
#' [pb_counts()] computed again: on an object saved
#' after a step was run again, run the later steps again.
#'
#' @return The updated `ScMaSigPro` object.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' path <- tempfile(fileext = ".rds")
#' saveRDS(scmp.ob, path)
#'
#' # Later, possibly with a newer version of scMaSigPro
#' ob <- updateObject(readRDS(path))
#' unlink(path)
#'
#' @seealso [ScMaSigPro-class]
#'
#' @importFrom S4Vectors updateObject
#' @importFrom methods .hasSlot
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @family classes
#'
#' @export
setMethod("updateObject", "ScMaSigPro", function(object, ...,
                                                 verbose = FALSE) {
    # Saved before the version slot existed (0.0.5 and earlier)
    if (!.hasSlot(object, "version")) {
        .check(isTRUE(validObject(object@Dense, test = TRUE)),
            msg = .rebin_message(object)
        )
        if (verbose) {
            message(
                "[updateObject] Adding the 'version' slot to an ScMaSigPro ",
                "object saved by an earlier version of scMaSigPro."
            )
        }
        object@version <- .scmp_version()
    }
    # Saved before the count filter of sc.p.vector()
    params <- object@Parameters
    missing_params <- c("min_counts", "min_bins_expressed")
    missing_params <- missing_params[
        !vapply(missing_params, .hasSlot, logical(1), object = params)
    ]
    if (length(missing_params) > 0L) {
        if (verbose) {
            message(
                "[updateObject] Adding the parameters ",
                paste(missing_params, collapse = " and "), " (0)."
            )
        }
        for (s in missing_params) {
            attr(params, s) <- 0
        }
        object@Parameters <- params
        object@version <- .scmp_version()
    }
    # Saved before the test of sc.p.vector() could be chosen: its fits used
    # the chi-squared test
    params <- object@Parameters
    if (!.hasSlot(params, "test")) {
        if (verbose) {
            message("[updateObject] Adding the parameter test (\"chisq\").")
        }
        attr(params, "test") <- "chisq"
        object@Parameters <- params
        object@version <- .scmp_version()
    }
    # Saved with the log_offset parameter, which chose the scale of the
    # offset: with log_offset = FALSE the fits used the size factors
    params <- object@Parameters
    if (.hasSlot(params, "log_offset")) {
        if (verbose) {
            message("[updateObject] Removing the parameter log_offset.")
        }
        fitted <- isTRUE(params@offset) &&
            isTRUE(any(object@Design@offset != 0))
        raw <- fitted && !isTRUE(attr(params, "log_offset"))
        attr(params, "log_offset") <- NULL
        object@Parameters <- params
        if (!fitted) {
            object@Design@offset <- numeric(0)
        } else if (raw) {
            object@Design@offset <- log(object@Design@offset)
            warning(.raw_offset_message(), call. = FALSE)
        }
        object@version <- .scmp_version()
    }
    # Fitted with a relabelled family
    if (.relabelled(object@Parameters@distribution)) {
        warning(.relabelled_message(object@Parameters@distribution),
            call. = FALSE
        )
    }
    object@Sparse <- updateObject(object@Sparse, ..., verbose = verbose)
    object@Dense <- updateObject(object@Dense, ..., verbose = verbose)
    validObject(object)
    # Saved after a step was run again by a version that kept the slots of
    # the later steps
    stale <- .stale_message(object)
    if (!is.null(stale)) {
        warning(stale, call. = FALSE)
    }
    object
})

# The warning of updateObject() for an object with a slot that was not
# computed from the slots before it, or NULL. The checks follow the
# workflow and the first that fails names the step to run again: the
# design has the binned pseudotime of the bins of the Dense slot (its first
# column; the bins may have been renamed or reordered since); tests made
# with an offset have it in the design; the models are those of non-flat
# genes, with a coefficient per column of the design; the genes of the
# lists have a model; the genes of the clusters are in the lists.
.stale_message <- function(object) {
    design <- object@Design@predictor_matrix
    other_bins <- FALSE
    if (nrow(design) > 0L) {
        alloc <- object@Design@assignment_matrix
        bin_time <- colData(object@Dense)[[colnames(alloc)[1L]]]
        other_bins <- !identical(
            sort(as.numeric(alloc[, 1L])), sort(as.numeric(bin_time))
        )
    }
    sol <- object@Estimate@significance_matrix
    coefs <- object@Estimate@coefficient_matrix
    genes <- unlist(object@Significant@genes, use.names = FALSE)
    clustered <- names(object@Significant@clusters)
    terms <- c("beta0", paste0("beta", colnames(design)))
    other_models <- !all(rownames(sol) %in% object@Profile@non_flat) ||
        !identical(colnames(coefs), terms)
    stale <- c(
        "sc.set.poly()" = other_bins,
        "sc.p.vector()" = length(object@Profile@p_values) > 0L &&
            isTRUE(object@Parameters@offset) &&
            length(object@Design@offset) == 0L,
        "sc.t.fit()" = nrow(coefs) > 0L && other_models,
        "sc.filter()" = !all(genes %in% rownames(sol)),
        "sc.cluster.trend()" = !all(clustered %in% genes)
    )
    if (!any(stale)) {
        return(NULL)
    }
    reasons <- c(
        "its design was built for other bins",
        "its tests used an offset that its design does not hold",
        "its models are those of other genes or of another design",
        "its gene lists hold genes without a model",
        "its clusters hold genes that are not in its gene lists"
    )
    first <- which(stale)[1L]
    paste0(
        "Some results of this object were not computed from its earlier ",
        "slots: ", reasons[first], ". Earlier versions of scMaSigPro kept ",
        "the results of the later steps when a step was run again: run ",
        names(stale)[first], " and the later steps on the object again."
    )
}

# The error of updateObject() for an object binned by scMaSigPro 0.0.4: the
# command that bins its cells again, with its own column names.
.rebin_message <- function(object) {
    paste0(
        "This ScMaSigPro object was binned by scMaSigPro 0.0.4: its 'Dense' ",
        "slot is not a valid SingleCellExperiment, so updateObject() cannot ",
        "update it. Bin the cells of its 'Sparse' slot again and re-run the ",
        "later steps; for this object as 'old_ob':\n",
        "new_ob <- sc.squeeze(as_scmp(old_ob@Sparse, ptime_col = \"",
        object@Parameters@ptime_col, "\", path_col = \"",
        object@Parameters@path_col, "\"))"
    )
}

# The warning of updateObject() for an object whose fits used the size
# factors as the offset (log_offset = FALSE)
.raw_offset_message <- function() {
    paste(
        "The fits of this object used the size factors of the bins as the",
        "offset (log_offset = FALSE). updateObject() stores their logarithm,",
        "the offset of sc.p.vector() and sc.t.fit(), but the 'Profile' and",
        "'Estimate' slots hold the fits made with the size factors: run",
        "sc.p.vector() and the later steps on the object again."
    )
}
