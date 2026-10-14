#' Polynomial design of the maSigPro model and the coefficients of a path
#'
#' Internal functions behind `sc.set.poly()`, `sc.t.fit()` and `plotTrend()`.
#' `.poly_design()` builds the regression design of the maSigPro model
#' (Conesa et al. 2006): every path is compared with the first one, the
#' reference, through a 0/1 dummy and through the products of that dummy
#' with each power of time up to `degree`. `.path_coeffs()` turns the
#' coefficients of such a model back into the polynomial of one path. They
#' were written from the published method and the maSigPro help pages of
#' `make.design.matrix()` and `reg.coeffs()`, then checked against those
#' functions as black boxes; no maSigPro code was used.
#'
#' `.poly_design()` takes the edesign that `sc.set.poly()` builds, a numeric
#' matrix: column 1 holds the binned pseudotime, column 2 the replicate
#' index, which the design does not use, and columns 3 onwards a 0/1
#' indicator for each path, the reference first. Every row belongs to one
#' path. With `<t>` the name of column 1, the columns of `dis` are the dummy
#' `<path>vs<ref>` of every other path, then for each power `k` of time the
#' term `<t>k` (`<t>` for `k = 1`) followed by `<t>kx<path>` for every other
#' path. In `groups.vector`, a power of time is labelled with the reference
#' and the other columns with `<path>vs<ref>`. With one path, `dis` holds the
#' powers of time only. The first power is the time column itself, so it
#' keeps the storage mode of `edesign`; higher powers are computed with `^`.
#' Unlike maSigPro, a design with a single row at the earliest time works.
#' (maSigPro also takes rows that belong to several groups, a common
#' starting point, and then leaves out dummies; `sc.set.poly()` never builds
#' such rows.)
#'
#' `.path_coeffs()` takes the coefficients of one gene, the intercept first
#' and then one per column of `dis`, and the `groups.vector` of the design,
#' with or without a label for the intercept. The reference is the shortest
#' label. The intercept and the coefficients labelled with the reference are
#' the reference path's polynomial, from the constant upwards; for another
#' path, the coefficients labelled `<path>vs<ref>` are added term by term. A
#' path with no dummy of its own starts from the reference intercept, and a
#' coefficient that could not be estimated (`NA`) stays `NA`.
#'
#' @param edesign Numeric matrix: time, replicate and one indicator column
#'   per path, with the rows and path columns named.
#' @param degree Degree of the polynomial, at least 1; its whole part is used.
#' @param coefficients Coefficients of one gene: a numeric vector or a
#'   one-row data frame, the intercept first.
#' @param groups.vector Path label of every column of `dis`, optionally
#'   preceded by the label of the intercept.
#' @param group Name of the path.
#'
#' @return `.poly_design()`: a list with the design matrix `dis`, its
#'   `groups.vector` and `edesign` unchanged. `.path_coeffs()`: an unnamed
#'   numeric vector, the coefficient of each power of time from 0 to
#'   `degree`.
#'
#' @references Conesa, A., Nueda, M.J., Ferrer, A. and Talon, M. (2006).
#' Bioinformatics 22, 1096-1102. \doi{10.1093/bioinformatics/btl056}
#'
#' @name poly_design
#' @noRd
NULL

.poly_design <- function(edesign, degree) {
    .check(
        is.numeric(degree) && length(degree) == 1L && is.finite(degree) &&
            degree >= 1,
        msg = "'poly_degree' must be a number of at least 1"
    )
    # The whole part, as sc.set.poly() records it
    degree <- as.integer(degree)
    paths <- colnames(edesign)[-c(1L, 2L)]
    ref <- paths[1L]
    others <- paths[-1L]
    # sprintf(), not paste0(): with one path `others` is empty and so is `vs`
    vs <- sprintf("%svs%s", others, ref)
    time <- edesign[, 1L]
    time_name <- colnames(edesign)[1L]
    # Indicator columns of the other paths, by position (columns 4 onwards)
    dummies <- lapply(seq_along(others) + 3L, function(j) edesign[, j])
    cols <- dummies
    col_names <- vs
    for (k in seq_len(degree)) {
        x <- if (k == 1L) time else time^k
        term <- if (k == 1L) time_name else paste0(time_name, k)
        cols <- c(cols, list(x), lapply(dummies, function(d) x * d))
        col_names <- c(col_names, term, sprintf("%sx%s", term, others))
    }
    dis <- matrix(unlist(cols, use.names = FALSE),
        nrow = nrow(edesign),
        dimnames = list(rownames(edesign), col_names)
    )
    list(
        dis = dis,
        groups.vector = c(vs, rep(c(ref, vs), degree)),
        edesign = edesign
    )
}

.path_coeffs <- function(coefficients, groups.vector, group) {
    coefficients <- as.numeric(coefficients)
    ref <- groups.vector[nchar(groups.vector) == min(nchar(groups.vector))][1L]
    # One label per coefficient after the intercept
    labels <- groups.vector
    if (length(labels) == length(coefficients)) {
        labels <- labels[-1L]
    }
    terms <- coefficients[-1L]
    base <- c(coefficients[1L], terms[labels == ref])
    shift <- terms[labels == paste0(group, "vs", ref)]
    if (length(shift) < length(base)) {
        shift <- c(0, shift)
    }
    base + shift
}
