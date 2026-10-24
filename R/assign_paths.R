#' @title Assign cells to paths from lineage weights
#'
#' @description
#' `assign_paths()` turns the lineage weights and pseudotimes of a trajectory
#' tool, one column per lineage as slingshot and tradeSeq give them, into
#' what scMaSigPro needs: one pseudotime value and one path label per cell.
#' Cells on every lineage are labelled as the root, and cells on some of the
#' lineages get the names of those lineages joined by `shared_sep`: the root
#' and link paths that [sc.restruct()] takes.
#'
#' @param pseudotime The pseudotime of each cell: a numeric matrix with the
#' row and column names of `weights`, such as `slingshot::slingPseudotime()`;
#' or, for a pseudotime that every lineage shares, as Monocle3 gives it, a
#' numeric vector with one value per row of `weights`, in that order or
#' named by cell.
#' @param weights A matrix or data frame of non-negative weights with one row
#' per cell, named by cell, and one column per lineage, named by lineage and
#' at least two of them, such as `slingshot::slingCurveWeights()`, the cell
#' weights of tradeSeq or the result of [paths_from_graph()].
#' @param min_weight A single positive number: a cell is on a lineage when its
#' weight there is at least `min_weight` and its pseudotime there is a finite
#' number. (Default is 0.5)
#' @param root_label The path label of the cells on every lineage.
#' (Default is "root")
#' @param shared_sep The separator of the lineage names in the label of the
#' cells on some but not all lineages; pass it to [sc.restruct()] as
#' `link_sep`. (Default is "_links_", the default of [sc.restruct()])
#' @param verbose Print a message when cells are left out. (Default is TRUE)
#'
#' @details
#' Each cell is placed on the lineages where its weight is at least
#' `min_weight` and its pseudotime is known. A cell on one lineage gets the
#' name of that lineage as its path, a cell on every lineage gets
#' `root_label`, and a cell on two or more lineages but not all of them gets
#' the names of those lineages, in the order of the columns of `weights`,
#' joined by `shared_sep`. A cell on no lineage is left out. The pseudotime
#' of a cell is the mean of its pseudotimes on its lineages, or its value
#' when `pseudotime` is a vector.
#'
#' The weights are compared with `min_weight` as they are given. Curve
#' weights of slingshot are close to 1 on every lineage for the cells of the
#' root; for weights that sum to 1 in each cell, as with
#' `slingCurveWeights(as.probs = TRUE)`, use a `min_weight` below
#' 1 / (number of lineages) or the root cells are left out. Soft weights can
#' place a few cells on a set of lineages that share no segment of the
#' trajectory beyond the root, such as two of three lineages that split at
#' the root; their label is built the same way, so relabel or remove them
#' before [sc.restruct()] if they are not a link.
#'
#' @section Input contract:
#' scMaSigPro needs one pseudotime value and one path label per cell, and
#' nothing else from the tool that inferred the trajectory. [create_scmp()]
#' takes them as a data frame named by cell, such as the result of
#' `assign_paths()`, and [as_scmp()] as two `colData()` columns of a
#' `SingleCellExperiment`. Code outside this package that reads a
#' trajectory, such as a companion package for Monocle3 objects or an
#' interactive branch selector, returns that data frame, or the R code that
#' builds it with [paths_from_graph()] and `assign_paths()`. Its path labels
#' follow the rules in Details, so that [sc.restruct()] can rebuild the
#' paths that share a root.
#'
#' @return A data frame with one row per cell that is on at least one
#' lineage, in the order of the rows of `weights` and named by cell, and two
#' columns: `Pseudotime` (numeric) and `Path` (character). Pass it to
#' [create_scmp()] as `cell_data` with `ptime_col = "Pseudotime"` and
#' `path_col = "Path"`, or add its columns to the `colData()` of a
#' `SingleCellExperiment` with the same cells for [as_scmp()], whose default
#' column names they are.
#'
#' @seealso [paths_from_graph()] for weights from a trajectory graph;
#' [sc.restruct()] for the root and link paths.
#'
#' @examples
#' # Three lineages: cells 1 and 2 are on all of them, cell 3 on A and B
#' weights <- matrix(
#'     c(
#'         1, 1, 1, 1, 0, 0,
#'         1, 1, 1, 0, 1, 0,
#'         1, 0.9, 0.2, 0, 0, 1
#'     ),
#'     ncol = 3,
#'     dimnames = list(paste0("cell", 1:6), c("A", "B", "C"))
#' )
#' pseudotime <- c(0.1, 0.5, 1.2, 2.1, 2.3, 1.8)
#' assign_paths(pseudotime, weights)
#'
#' # A pseudotime for each lineage, as slingshot gives it
#' pseudotime <- cbind(A = pseudotime, B = pseudotime + 0.1, C = pseudotime)
#' rownames(pseudotime) <- rownames(weights)
#' assign_paths(pseudotime, weights, shared_sep = "|")
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
assign_paths <- function(pseudotime, weights, min_weight = 0.5,
                         root_label = "root", shared_sep = "_links_",
                         verbose = TRUE) {
    # The weights: cells as rows, at least two lineages as columns, all named
    .check(is.matrix(weights) || is.data.frame(weights),
        msg = paste(
            "'weights' must be a matrix or a data frame with one row per",
            "cell and one column per lineage."
        )
    )
    weights <- as.matrix(weights)
    .check(is.numeric(weights) && all(is.finite(weights)) && all(weights >= 0),
        msg = "'weights' must hold non-negative numbers, without NA."
    )
    .check(nrow(weights) >= 1L && ncol(weights) >= 2L,
        msg = "'weights' must have one row per cell and two or more lineages."
    )
    cells <- rownames(weights)
    lineages <- colnames(weights)
    .check(.unique_names(cells),
        msg = "The rows of 'weights' must be named by cell, each name once."
    )
    .check(.unique_names(lineages),
        msg = paste(
            "The columns of 'weights' must be named by lineage, each name",
            "once."
        )
    )
    .check(
        is.numeric(min_weight) && length(min_weight) == 1L &&
            is.finite(min_weight) && min_weight > 0,
        msg = "'min_weight' must be a single positive number."
    )
    .check_strings(root_label = root_label, shared_sep = shared_sep)
    .check(nzchar(root_label) && nzchar(shared_sep),
        msg = "'root_label' and 'shared_sep' must not be empty strings."
    )
    # Every label must say which lineages a cell is on
    .check(
        !root_label %in% lineages &&
            !any(grepl(shared_sep, c(root_label, lineages), fixed = TRUE)),
        msg = paste(
            "The lineage names must differ from 'root_label', and neither",
            "may contain 'shared_sep'."
        )
    )

    # The pseudotime: a matrix like the weights, or one value per cell
    if (is.data.frame(pseudotime)) {
        pseudotime <- as.matrix(pseudotime)
    }
    shared_time <- !is.matrix(pseudotime)
    if (shared_time) {
        time_names <- names(pseudotime)
        .check(
            is.numeric(pseudotime) && length(pseudotime) == length(cells) &&
                (is.null(time_names) || identical(time_names, cells)),
            msg = paste(
                "A 'pseudotime' vector must hold one number per row of",
                "'weights', unnamed or named by the same cells in the same",
                "order."
            )
        )
        time <- matrix(pseudotime, nrow = length(cells), ncol = ncol(weights))
    } else {
        .check(
            is.numeric(pseudotime) &&
                identical(rownames(pseudotime), cells) &&
                identical(colnames(pseudotime), lineages),
            msg = paste(
                "A 'pseudotime' matrix must be numeric, with the row and",
                "column names of 'weights' in the same order."
            )
        )
        time <- pseudotime
    }

    # The lineages of each cell
    on <- weights >= min_weight & is.finite(time)
    n_on <- rowSums(on)
    keep <- n_on > 0
    .check(any(keep),
        msg = paste(
            "No cell is on a lineage: every weight is below 'min_weight' or",
            "has no pseudotime."
        )
    )
    if (verbose && !all(keep)) {
        message(
            sum(!keep), " of ", length(cells), " cells are on no lineage ",
            "(no weight of at least ", min_weight, " with a pseudotime) and ",
            "are left out."
        )
    }
    on <- on[keep, , drop = FALSE]
    time <- time[keep, , drop = FALSE]

    path <- apply(on, 1L, function(x) paste(lineages[x], collapse = shared_sep))
    path[n_on[keep] == length(lineages)] <- root_label
    if (shared_time) {
        ptime <- pseudotime[keep]
    } else {
        ptime <- vapply(seq_len(nrow(on)), function(i) {
            mean(time[i, on[i, ]])
        }, numeric(1))
    }

    data.frame(
        Pseudotime = unname(ptime),
        Path = unname(path),
        row.names = cells[keep]
    )
}

# TRUE when x is a character vector of names, none missing or empty, each
# given once.
.unique_names <- function(x) {
    is.character(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
}
