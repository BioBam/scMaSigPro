# assign_paths(): lineage weights and pseudotimes to one pseudotime and one
# path label per cell, with known answers on small matrices.

# Three lineages: cells 1 and 2 on all of them, cell 3 on A and B, cells 4
# to 6 on one lineage each
ap_weights <- function() {
    matrix(
        c(
            1, 1, 1, 1, 0, 0,
            1, 1, 1, 0, 1, 0,
            1, 0.9, 0.2, 0, 0, 1
        ),
        ncol = 3,
        dimnames = list(paste0("cell", 1:6), c("A", "B", "C"))
    )
}
ap_time <- c(0.1, 0.5, 1.2, 2.1, 2.3, 1.8)

test_that("assign_paths() labels root, link and leaf cells", {
    out <- assign_paths(ap_time, ap_weights())
    expect_identical(out, data.frame(
        Pseudotime = ap_time,
        Path = c("root", "root", "A_links_B", "A", "B", "C"),
        row.names = paste0("cell", 1:6)
    ))
    # The same with weights as a data frame and a pseudotime named by cell
    named <- stats::setNames(ap_time, paste0("cell", 1:6))
    expect_identical(
        assign_paths(named, as.data.frame(ap_weights())), out
    )
})

test_that("min_weight, root_label and shared_sep set the labels", {
    # Cell 2 has weight 0.9 on C: not on C at 0.95
    out <- assign_paths(ap_time, ap_weights(),
        min_weight = 0.95,
        root_label = "Root", shared_sep = "|"
    )
    expect_identical(
        out$Path, c("Root", "A|B", "A|B", "A", "B", "C")
    )
    # Every weight counts at 0.1, so cell 3 is on all three lineages
    out <- assign_paths(ap_time, ap_weights(), min_weight = 0.1)
    expect_identical(out$Path[1:3], rep("root", 3))
})

test_that("two lineages give root and leaf labels only", {
    w <- ap_weights()[, c("A", "C")]
    out <- assign_paths(ap_time, w, verbose = FALSE)
    expect_identical(rownames(out), paste0("cell", c(1:4, 6)))
    expect_identical(out$Path, c("root", "root", "A", "A", "C"))
})

test_that("a pseudotime matrix gives the mean on the lineages of a cell", {
    w <- ap_weights()
    time <- cbind(A = ap_time, B = ap_time + 0.3, C = ap_time + 0.6)
    rownames(time) <- rownames(w)
    # Pseudotimes off the lineages of a cell may be NA, as in slingshot
    time[w == 0] <- NA
    out <- assign_paths(time, w)
    expect_equal(out$Pseudotime, c(
        mean(c(0.1, 0.4, 0.7)), mean(c(0.5, 0.8, 1.1)), mean(c(1.2, 1.5)),
        2.1, 2.6, 2.4
    ))
    expect_identical(out$Path, c("root", "root", "A_links_B", "A", "B", "C"))
    # Without a pseudotime on A, cell 4 is on no lineage
    time["cell4", "A"] <- NA
    expect_message(
        out <- assign_paths(time, w),
        paste(
            "1 of 6 cells are on no lineage (no weight of at least 0.5 with",
            "a pseudotime) and are left out."
        ),
        fixed = TRUE
    )
    expect_identical(rownames(out), paste0("cell", c(1:3, 5:6)))
    # Monocle3 gives Inf to cells it cannot reach from the root
    time["cell4", "A"] <- Inf
    expect_no_message(out <- assign_paths(time, w, verbose = FALSE))
    expect_false("cell4" %in% rownames(out))
})

test_that("cells on no lineage are left out with a message", {
    w <- rbind(ap_weights(), cell7 = c(0.2, 0.3, 0.1))
    expect_message(
        out <- assign_paths(c(ap_time, 3), w),
        "1 of 7 cells are on no lineage",
        fixed = TRUE
    )
    expect_identical(out, assign_paths(ap_time, ap_weights()))
    expect_no_message(assign_paths(c(ap_time, 3), w, verbose = FALSE))
})

test_that("assign_paths() checks its arguments", {
    w <- ap_weights()
    msgs <- list(
        weights = paste(
            "'weights' must be a matrix or a data frame with one row per",
            "cell and one column per lineage."
        ),
        values = "'weights' must hold non-negative numbers, without NA.",
        shape = paste(
            "'weights' must have one row per cell and two or more",
            "lineages."
        ),
        cells = "The rows of 'weights' must be named by cell, each name once.",
        lineages = paste(
            "The columns of 'weights' must be named by lineage, each name",
            "once."
        ),
        min_weight = "'min_weight' must be a single positive number.",
        labels = paste(
            "The lineage names must differ from 'root_label', and neither",
            "may contain 'shared_sep'."
        ),
        vector = paste(
            "A 'pseudotime' vector must hold one number per row of",
            "'weights', unnamed or named by the same cells in the same",
            "order."
        ),
        matrix = paste(
            "A 'pseudotime' matrix must be numeric, with the row and",
            "column names of 'weights' in the same order."
        ),
        none = paste(
            "No cell is on a lineage: every weight is below 'min_weight' or",
            "has no pseudotime."
        )
    )
    expect_error(assign_paths(ap_time, c(1, 0)), msgs$weights, fixed = TRUE)
    bad <- w
    bad[1, 1] <- NA
    expect_error(assign_paths(ap_time, bad), msgs$values, fixed = TRUE)
    bad[1, 1] <- -1
    expect_error(assign_paths(ap_time, bad), msgs$values, fixed = TRUE)
    expect_error(
        assign_paths(ap_time, w > 0), msgs$values,
        fixed = TRUE
    )
    expect_error(
        assign_paths(ap_time, w[, "A", drop = FALSE]), msgs$shape,
        fixed = TRUE
    )
    bad <- w
    rownames(bad) <- NULL
    expect_error(assign_paths(ap_time, bad), msgs$cells, fixed = TRUE)
    rownames(bad) <- c(paste0("cell", 1:5), "cell1")
    expect_error(assign_paths(ap_time, bad), msgs$cells, fixed = TRUE)
    bad <- w
    colnames(bad) <- c("A", "B", "A")
    expect_error(assign_paths(ap_time, bad), msgs$lineages, fixed = TRUE)
    colnames(bad) <- NULL
    expect_error(assign_paths(ap_time, bad), msgs$lineages, fixed = TRUE)
    for (mw in list(0, -1, NA_real_, c(0.5, 0.6), "0.5", Inf)) {
        expect_error(assign_paths(ap_time, w, min_weight = mw),
            msgs$min_weight,
            fixed = TRUE
        )
    }
    expect_error(assign_paths(ap_time, w, root_label = NA_character_),
        "'root_label' must be a single character string.",
        fixed = TRUE
    )
    expect_error(assign_paths(ap_time, w, shared_sep = c("|", "_")),
        "'shared_sep' must be a single character string.",
        fixed = TRUE
    )
    expect_error(assign_paths(ap_time, w, shared_sep = ""),
        "'root_label' and 'shared_sep' must not be empty strings.",
        fixed = TRUE
    )
    expect_error(assign_paths(ap_time, w, root_label = "A"),
        msgs$labels,
        fixed = TRUE
    )
    bad <- w
    colnames(bad) <- c("A", "B_links_C", "C")
    expect_error(assign_paths(ap_time, bad), msgs$labels, fixed = TRUE)
    expect_error(assign_paths(ap_time, w, root_label = "r|t", shared_sep = "|"),
        msgs$labels,
        fixed = TRUE
    )
    expect_error(assign_paths(ap_time[-1], w), msgs$vector, fixed = TRUE)
    expect_error(assign_paths(as.character(ap_time), w), msgs$vector,
        fixed = TRUE
    )
    expect_error(
        assign_paths(stats::setNames(ap_time, rev(rownames(w))), w),
        msgs$vector,
        fixed = TRUE
    )
    time <- w
    expect_error(assign_paths(time[, 3:1], w), msgs$matrix, fixed = TRUE)
    expect_error(assign_paths(unname(time), w), msgs$matrix, fixed = TRUE)
    expect_error(assign_paths(ap_time, w, min_weight = 2), msgs$none,
        fixed = TRUE
    )
})

test_that("the labels of assign_paths() run through sc.restruct()", {
    # Synthetic counts: 30 root cells, 20 cells on the A-B link, 25 cells on
    # each of A, B and C, 20 genes
    set.seed(11)
    n <- c(root = 30, link = 20, A = 25, B = 25, C = 25)
    cells <- paste0("c", seq_len(sum(n)))
    group <- rep(names(n), n)
    w <- cbind(
        A = as.numeric(group %in% c("root", "link", "A")),
        B = as.numeric(group %in% c("root", "link", "B")),
        C = as.numeric(group %in% c("root", "C"))
    )
    rownames(w) <- cells
    time <- ifelse(group == "root", 0, ifelse(group == "link", 1, 2)) +
        stats::runif(length(cells))
    cell_data <- assign_paths(time, w)
    expect_identical(
        as.vector(table(cell_data$Path)[c("root", "A_links_B", "A", "B", "C")]),
        as.integer(n)
    )
    counts <- matrix(stats::rpois(20 * length(cells), 5),
        nrow = 20,
        dimnames = list(paste0("g", 1:20), cells)
    )
    ob <- create_scmp(counts, cell_data,
        ptime_col = "Pseudotime", path_col = "Path"
    )
    ob <- suppressMessages(sc.squeeze(ob, verbose = FALSE))
    ob <- sc.restruct(ob,
        end_node_list = list("A", "B", "C"), root_node = "root",
        link_node_list = list("A_links_B"), verbose = FALSE
    )
    expect_valid_scmp(ob)
    expect_setequal(cDense(ob)$Path, c("A", "B", "C"))
})
