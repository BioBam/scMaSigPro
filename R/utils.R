#' @title Gene lists as a 0/1 table
#'
#' @description
#' `fromListWithNames()` turns a named list of gene lists into a data frame
#' with one row per gene, named by gene, and a 0/1 column per list, as
#' `UpSetR::fromList()` does, with the row names.
#'
#' @param input A named list of character vectors.
#'
#' @return A data frame of 0 and 1, one row per gene and one column per list.
#'
#' @keywords internal
#'
#' @seealso [plotIntersect()], which calls it.
fromListWithNames <- function(input) {
    # Get unique elements (genes)
    elements <- unique(unlist(input))

    # Create binary matrix
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = FALSE))

    # Filter out rows with no data
    data <- data[which(rowSums(data) != 0), ]

    # Set column names
    names(data) <- names(input)

    # Set row names
    rownames(data) <- elements[which(rowSums(data) != 0)]

    return(data)
}

#' @title Colours of the paths
#'
#' @description
#' `scmp_colors()` returns `n` contrasting colours: those of a fixed palette
#' of nine, then named colours of [grDevices::colors()].
#'
#' @param n The number of colours.
#'
#' @return A character vector of `n` colours.
#'
#' @keywords internal
#'
#' @seealso [plotTrend()] and [plotTrendCluster()], which call it.
scmp_colors <- function(n) {
    # Define Palette
    ten_pal <- c("#d95f02", "#4daf4a", "#377eb8", "#E69F00", "#f781bf", "#56B4E9", "#a65628", "#009E73", "#7570b3")
    extra_pal <- unique(grDevices::colors(distinct = TRUE))
    extra_pal <- extra_pal[grep("^[a-zA-Z]+$", extra_pal)]
    extra_pal <- extra_pal[c(2:length(extra_pal))]
    extra_pal <- extra_pal[!extra_pal %in% ten_pal]

    if (n <= 9) {
        return(ten_pal[1:n])
    } else {
        return(c(ten_pal, extra_pal[1:(n - 9)]))
    }
}

#' @title Summarise a Column Within Groups
#'
#' @description
#' `.summarise_by()` applies `fun` to column `value` of `df` within each group
#' of rows that share the values of the columns `by`, as
#' `dplyr::group_by()` followed by `dplyr::summarise(.groups = "drop")` did:
#' one row per group, ordered by the `by` columns in turn (radix order, so
#' numbers increase and strings sort in the C locale, as dplyr sorts groups),
#' with the columns `by` and then `value`. Keys are compared exactly, `NA` is a
#' group of its own (sorted last), and the values reach `fun` in row order.
#'
#' @param df A data frame.
#' @param by Names of the grouping columns.
#' @param value Name of the column to summarise.
#' @param fun A function returning one number.
#' @param ... Further arguments to `fun`.
#'
#' @return A data frame with automatic row names.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @noRd
.summarise_by <- function(df, by, value, fun, ...) {
    keys <- unname(as.list(df[by]))
    sorted <- df[do.call(order, c(keys, method = "radix")), , drop = FALSE]
    n <- nrow(sorted)

    # A group starts at the first row and wherever a key changes
    starts <- rep(TRUE, n)
    if (n > 1) {
        same <- lapply(sorted[by], function(k) {
            a <- k[-1]
            b <- k[-n]
            (!is.na(a) & !is.na(b) & a == b) | (is.na(a) & is.na(b))
        })
        starts[-1] <- !Reduce(`&`, same)
    }

    out <- sorted[starts, by, drop = FALSE]
    out[[value]] <- vapply(split(sorted[[value]], cumsum(starts)), fun,
        numeric(1), ...,
        USE.NAMES = FALSE
    )
    rownames(out) <- NULL
    out
}
