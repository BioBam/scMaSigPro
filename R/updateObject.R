#' @title Update an ScMaSigPro object to the current class definition
#'
#' @description
#' `updateObject()` brings an \code{\link{ScMaSigPro}} object saved by an
#' earlier version of \pkg{scMaSigPro}, for example with `saveRDS()` or
#' `save()`, up to date. Call it on such an object after `readRDS()` or
#' `load()`, before the next step of the workflow.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param ... Additional arguments, passed to the `updateObject()` methods of
#' the `Sparse` and `Dense` slots.
#' @param verbose Print a message when the object is updated. (Default is
#' FALSE)
#'
#' @details
#' Objects saved before the `version` slot was added (scMaSigPro 0.0.5 and
#' earlier) fail `validObject()`; `updateObject()` adds the slot with the
#' version of the package that is running. For every object, the `Sparse`
#' and `Dense` slots are passed to the `updateObject()` method of
#' \pkg{SingleCellExperiment}, and the result is checked with
#' `validObject()`. A current object, such as one created by this version
#' of the package, is returned unchanged.
#'
#' Objects saved before the `min_counts` and `min_bins_expressed`
#' parameters existed get both set to 0, the values that keep every gene
#' ([sc.p.vector()] had no such filter then), and the version of the
#' package that is running.
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
#' @seealso \code{\link{ScMaSigPro}}
#'
#' @importFrom S4Vectors updateObject
#' @importFrom methods .hasSlot
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
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
    # Saved before the count filter of sc.p.vector() (decision D17)
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
    # Fitted with a relabelled family (ledger LK-4)
    if (.relabelled(object@Parameters@distribution)) {
        warning(.relabelled_message(object@Parameters@distribution),
            call. = FALSE
        )
    }
    object@Sparse <- updateObject(object@Sparse, ..., verbose = verbose)
    object@Dense <- updateObject(object@Dense, ..., verbose = verbose)
    validObject(object)
    object
})

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
