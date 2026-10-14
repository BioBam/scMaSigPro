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
