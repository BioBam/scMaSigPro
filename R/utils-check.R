#' Stop unless a condition is TRUE
#'
#' Internal replacement for `assertthat::assert_that(cond, msg = msg)`. The
#' check passes only when `cond` is a single `TRUE` ([isTRUE()]); `FALSE`,
#' `NA`, a non-logical value and a vector of any other length fail. The
#' error has no call, so it prints as `Error: <msg>`, as assert_that() did.
#' `msg` is evaluated only when the check fails; several strings are joined
#' with newlines.
#'
#' @param cond The condition to check.
#' @param msg A character vector: the error message.
#'
#' @return `invisible(TRUE)` if the check passes; otherwise it signals an
#'   error.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @noRd
.check <- function(cond, msg) {
    if (!isTRUE(cond)) {
        stop(paste(msg, collapse = "\n"), call. = FALSE)
    }
    invisible(TRUE)
}

#' Stop unless an object can go through a function of the package
#'
#' `.check_scmp()` opens every exported function that takes an
#' `ScMaSigPro` object, before any work is done: the object must be of that
#' class and have the `version` slot, which objects saved by scMaSigPro
#' 0.0.5 and earlier lack until [updateObject()] adds it, and the
#' `min_counts` and `min_bins_expressed` parameters, which [updateObject()]
#' adds to objects saved before them, but not the `log_offset` parameter,
#' which [updateObject()] removes.
#' `.check_dense_counts()` opens the functions that read the pseudo-bulk
#' counts: the `Dense` slot must hold a `bulk.counts` assay, which it lacks
#' after [cDense<-] changed the number of bins and before [pb_counts()] was
#' run. Given `bins`, the bins of the design, the assay must also have a
#' column for each of them, as the fits select the columns by name.
#'
#' @param scmpObj The object given to the function.
#' @param bins `NULL`, or the names of the bins the function reads.
#'
#' @return `invisible(TRUE)` if the check passes; otherwise it signals an
#'   error.
#'
#' @noRd
.check_scmp <- function(scmpObj) {
    .check(is(scmpObj, "ScMaSigPro"),
        msg = "Please provide object of class 'ScMaSigPro'"
    )
    .check(.scmp_current(scmpObj),
        msg = paste(
            "This ScMaSigPro object was saved by an earlier version of",
            "scMaSigPro: update it with updateObject() first."
        )
    )
}

# TRUE when the object has every slot of the current class definition that
# updateObject() adds to objects saved by earlier versions, and none of the
# slots it removes
.scmp_current <- function(object) {
    .hasSlot(object, "version") &&
        .hasSlot(object@Parameters, "min_counts") &&
        .hasSlot(object@Parameters, "min_bins_expressed") &&
        !.hasSlot(object@Parameters, "log_offset")
}

.check_dense_counts <- function(scmpObj, bins = NULL) {
    dense <- scmpObj@Dense
    .check("bulk.counts" %in% assayNames(dense),
        msg = paste(
            "The 'Dense' slot has no 'bulk.counts' assay: run pb_counts()",
            "after changing the bins with cDense<-."
        )
    )
    if (!is.null(bins)) {
        counts <- assay(dense, "bulk.counts", withDimnames = FALSE)
        .check(all(bins %in% colnames(counts)),
            msg = paste(
                "The 'bulk.counts' assay has no counts for some bins of the",
                "design: run pb_counts() after changing the bins with",
                "cDense<-, then sc.set.poly()."
            )
        )
    }
}

#' Check the arguments that take a choice or a name
#'
#' `.check_choice()` returns `x` if it is one of `choices`, given once, and
#' otherwise stops with `msg`, the message of the calling function. Unlike
#' [match.arg()], it does not match partially and keeps that message.
#'
#' `.check_strings()` stops unless each of its named arguments is a single
#' non-missing character string, such as a column name or a label.
#'
#' @param x The value of the argument.
#' @param choices A character vector, the accepted values.
#' @param msg A character vector: the error message.
#' @param ... Named arguments to check.
#'
#' @return `.check_choice()`: `x`. `.check_strings()`: `invisible(TRUE)`.
#'   Both signal an error if a check fails.
#'
#' @noRd
.check_choice <- function(x, choices, msg) {
    .check(length(x) == 1L && x %in% choices, msg = msg)
    x
}

.check_strings <- function(...) {
    args <- list(...)
    for (name in names(args)) {
        x <- args[[name]]
        .check(is.character(x) && length(x) == 1L && !is.na(x),
            msg = paste0("'", name, "' must be a single character string.")
        )
    }
    invisible(TRUE)
}
