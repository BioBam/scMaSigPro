#' Stepwise selection of the terms of a gene's regression model
#'
#' Internal functions behind `sc.t.fit(selection_method = ...)`. Each selects
#' the significant terms of the GLM `y ~ .` on the columns of `d` with one of
#' the stepwise strategies of maSigPro (Conesa et al. 2006; Nueda et al.
#' 2014) and returns the final [glm()] fit. They were written from the
#' published method and the maSigPro help pages, then checked against the
#' maSigPro functions as black boxes; no maSigPro code was used.
#'
#' Terms are handled by column position, so names with spaces or `|` need no
#' cleaning. A term's p-value is that of its coefficient in `summary()`
#' (column 4); aliased terms have none. Every fit passes `family`,
#' `useWeight`, `useOffset`, `epsilon` and `max_it` to [glm()].
#'
#' * `sc.stepback()`: remove the term with the largest p-value while it is
#'   above `alfa`; aliased terms stay. Removing the last estimable term gives
#'   the intercept-only model.
#' * `sc.stepfor()`: the candidates are the terms estimable in the full
#'   model; the one with the smallest p-value when added enters if that
#'   p-value is below `alfa`.
#' * `sc.two.ways.stepback()`: as `sc.stepback()`, but after each removal
#'   the removed term with the smallest p-value when added back re-enters if
#'   that p-value is at most `alfa`, one at a time. Aliased terms are dropped.
#'   When a state comes back (a cycle), the re-entries stop and the worst
#'   term is removed while above `alfa`, as in `sc.stepback()`, so every term
#'   of the returned model is significant.
#' * `sc.two.ways.stepfor()`: as `sc.stepfor()` with "at most `alfa`"; after
#'   each entry the worst term is removed while above `alfa` and queued as a
#'   candidate again. The fit in which the last candidate enters is returned
#'   unchecked; removing every term gives the intercept-only model.
#'
#' Ties go to the first term. With no p-value at all, [max()] and [min()]
#' warn as in maSigPro. Where maSigPro loops for ever, the two-ways methods
#' stop when a state comes back: forward with the current model, backward
#' after the final removals described above. The forward methods keep
#' the row names of `d`, so the influential observations of their fits are
#' named by bin (maSigPro numbers the rows, so its influence matrix lists the
#' genes but every entry is NA); as in maSigPro, the backward methods keep the
#' row names of `d` until one column is left.
#'
#' @param y Numeric response, one value per row of `d`.
#' @param d A data frame of numeric predictors, one column per term.
#' @param alfa Significance level for a term to enter or to stay.
#' @param family,epsilon,max_it Passed to [glm()] (`maxit = max_it`).
#' @param useOffset,useWeight `NULL` or one offset or weight per row of `d`.
#'
#' @return The final `glm` fit, `y ~ 1` when no term is selected.
#'
#' @references Conesa, A., Nueda, M.J., Ferrer, A. and Talon, M. (2006).
#' Bioinformatics 22, 1096-1102. \doi{10.1093/bioinformatics/btl056}
#'
#' Nueda, M.J., Tarazona, S. and Conesa, A. (2014). Bioinformatics 30,
#' 2598-2602. \doi{10.1093/bioinformatics/btu333}
#'
#' @importFrom stats gaussian glm
#' @name stepwise
#' @noRd
NULL

# Fit y ~ . on the columns `cols` (positions in d), or y ~ 1 when there are
# none. reset_rows numbers the rows instead of keeping the names of d.
.term_fit <- function(y, d, cols, family, epsilon, useOffset, useWeight,
                      max_it, reset_rows = FALSE) {
    if (length(cols) == 0L) {
        return(glm(y ~ 1,
            family = family, weights = useWeight, offset = useOffset,
            epsilon = epsilon, maxit = max_it
        ))
    }
    design <- d[, cols, drop = FALSE]
    if (reset_rows) {
        rownames(design) <- NULL
    }
    glm(y ~ .,
        family = family, data = design, weights = useWeight,
        offset = useOffset, epsilon = epsilon, maxit = max_it
    )
}

# TRUE for every term (column of the fit's data) whose coefficient is
# estimable; one coefficient per numeric column, after the intercept.
.term_estimable <- function(fit) !is.na(fit$coefficients[-1L])

# p-value of every term of the fit, in column order; NA when not estimable
.term_p <- function(fit) {
    est <- .term_estimable(fit)
    p <- rep(NA_real_, length(est))
    p[est] <- summary(fit)$coefficients[-1L, 4L]
    p
}

# p-value of the last term of the fit: the term just added
.last_p <- function(fit) {
    p <- .term_p(fit)
    p[length(p)]
}

# p-value of each candidate when it is added after the selected terms
.add_p <- function(fit_cols, sel, cand) {
    vapply(cand, function(j) .last_p(fit_cols(c(sel, j))), numeric(1))
}

sc.stepback <- function(y, d, alfa = 0.05, family = gaussian(),
                        epsilon = 1e-05, useOffset = NULL, useWeight = NULL,
                        max_it = 25) {
    fit_cols <- function(cols) {
        .term_fit(y, d, cols, family, epsilon, useOffset, useWeight, max_it,
            reset_rows = length(cols) == 1L && ncol(d) > 1L
        )
    }
    cols <- seq_len(ncol(d))
    repeat {
        fit <- fit_cols(cols)
        p <- .term_p(fit)
        if (max(p, na.rm = TRUE) <= alfa) {
            return(fit)
        }
        if (sum(.term_estimable(fit)) == 1L) {
            return(fit_cols(integer(0)))
        }
        cols <- cols[-which.max(p)]
    }
}

sc.stepfor <- function(y, d, alfa = 0.05, family = gaussian(),
                       epsilon = 1e-05, useOffset = NULL, useWeight = NULL,
                       max_it = 25) {
    fit_cols <- function(cols) {
        .term_fit(y, d, cols, family, epsilon, useOffset, useWeight, max_it,
            reset_rows = FALSE
        )
    }
    cand <- which(.term_estimable(fit_cols(seq_len(ncol(d)))))
    sel <- integer(0)
    while (length(cand) > 0L) {
        p <- .add_p(fit_cols, sel, cand)
        if (min(p, na.rm = TRUE) >= alfa) {
            break
        }
        best <- which.min(p)
        sel <- c(sel, cand[best])
        cand <- cand[-best]
    }
    fit_cols(sel)
}

sc.two.ways.stepback <- function(y, d, alfa = 0.05, family = gaussian(),
                                 epsilon = 1e-05, useOffset = NULL,
                                 useWeight = NULL, max_it = 25) {
    reset_rows <- FALSE
    fit_cols <- function(cols) {
        .term_fit(y, d, cols, family, epsilon, useOffset, useWeight, max_it,
            reset_rows = reset_rows
        )
    }
    cols <- seq_len(ncol(d))
    out <- integer(0)
    seen <- character(0)
    fit <- fit_cols(cols)
    repeat {
        p <- .term_p(fit)
        if (max(p, na.rm = TRUE) <= alfa) {
            return(fit)
        }
        # Removals and re-entries can cycle, and maSigPro then never returns:
        # when the terms in and out of the model come back, stop the
        # re-entries and remove the worst term while it is above alfa
        state <- paste(c(cols, 0L, out), collapse = " ")
        if (state %in% seen) {
            repeat {
                worst <- which.max(p)
                cols <- cols[.term_estimable(fit) & seq_along(cols) != worst]
                if (length(cols) == 0L) {
                    return(fit_cols(integer(0)))
                }
                reset_rows <- reset_rows || length(cols) == 1L
                fit <- fit_cols(cols)
                p <- .term_p(fit)
                if (max(p, na.rm = TRUE) <= alfa) {
                    return(fit)
                }
            }
        }
        seen <- c(seen, state)
        worst <- which.max(p)
        out <- c(out, cols[worst])
        cols <- cols[.term_estimable(fit) & seq_along(cols) != worst]
        if (length(cols) == 0L) {
            return(fit_cols(integer(0)))
        }
        reset_rows <- reset_rows || length(cols) == 1L
        # Let the removed terms back in, the most significant first
        while (length(out) > 0L) {
            p_out <- .add_p(fit_cols, cols, out)
            if (all(is.na(p_out)) || min(p_out, na.rm = TRUE) > alfa) {
                break
            }
            best <- which.min(p_out)
            cols <- c(cols, out[best])
            out <- out[-best]
        }
        fit <- fit_cols(cols)
    }
}

sc.two.ways.stepfor <- function(y, d, alfa = 0.05, family = gaussian(),
                                epsilon = 1e-05, useOffset = NULL,
                                useWeight = NULL, max_it = 25) {
    fit_cols <- function(cols) {
        .term_fit(y, d, cols, family, epsilon, useOffset, useWeight, max_it,
            reset_rows = FALSE
        )
    }
    cand <- which(.term_estimable(fit_cols(seq_len(ncol(d)))))
    sel <- integer(0)
    seen <- character(0)
    while (length(cand) > 0L) {
        # Entries and removals can cycle, and maSigPro then never returns:
        # stop when a selection and candidate list come back
        state <- paste(c(sel, 0L, cand), collapse = " ")
        if (state %in% seen) {
            break
        }
        seen <- c(seen, state)
        p <- .add_p(fit_cols, sel, cand)
        if (min(p, na.rm = TRUE) > alfa) {
            break
        }
        best <- which.min(p)
        sel <- c(sel, cand[best])
        cand <- cand[-best]
        fit <- fit_cols(sel)
        if (length(cand) == 0L) {
            return(fit)
        }
        # Remove the terms that the new one made non-significant
        repeat {
            p_sel <- .term_p(fit)
            if (max(p_sel, na.rm = TRUE) <= alfa) {
                break
            }
            worst <- which.max(p_sel)
            cand <- c(cand, sel[worst])
            sel <- sel[-worst]
            if (length(sel) == 0L) {
                # No term is left (maSigPro stops with an error here)
                return(fit_cols(integer(0)))
            }
            fit <- fit_cols(sel)
        }
    }
    fit_cols(sel)
}
