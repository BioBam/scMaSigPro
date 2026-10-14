#' Clean or Check Backticks in a String
#'
#' This function either checks for the presence of backticks in a string or removes them, based on the specified action.
#'
#' @param input_string A character string to be processed.
#' @param action A character string specifying the action to perform. Either "check" to check for backticks or "remove" to remove backticks. Defaults to "check".
#'
#' @return If action is "check", returns a logical value indicating the presence of backticks. If action is "remove", returns the input string with all backticks removed.
#'
#' @keywords internal
clean_string <- function(input_string, action = c("check", "remove")) {
    action <- match.arg(action)

    if (action == "check") {
        return(grepl("`", input_string, fixed = TRUE))
    } else if (action == "remove") {
        # Remove Backticks
        processed_string <- gsub("`", "", input_string, fixed = TRUE)
        return(processed_string)
    }
}


#' Convert List of Named Vectors to Data Frame with Row Names for UpSetR
#'
#' This function converts a list of named vectors to a data frame compatible with UpSetR,
#' retaining the row names corresponding to the unique elements (genes) present in the input list.
#'
#' @param input A list of named vectors to be converted to a data frame compatible with UpSetR.
#'
#' @return A data frame with binary values indicating the presence (1) or absence (0) of each element
#' in the sets, and row names corresponding to the unique elements.
#'
#' @keywords internal
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

#' Get Color Palette
#'
#' This function returns a specified number of contrasting colors from a predefined palette.
#' If the requested number of colors exceeds the length of the predefined palette,
#' additional unique colors from the default R color set are included.
#'
#' @param n An integer specifying the number of colors to return.
#'
#' @return A character vector of hex color codes.
#'
#' @keywords internal
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
