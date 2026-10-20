# paths_from_graph(): a cell-to-node map and an edge list to 0/1 lineage
# weights, with known answers on a 10-node tree:
#
#   R - A - B - C - D
#        \   \
#         \   E - F
#          G - H - I
#
# The end nodes are D, F and I; B is shared by the paths to D and F, and R
# and A by all three.

pg_edges <- function() {
    rbind(
        c("R", "A"), c("A", "B"), c("B", "C"), c("C", "D"), c("B", "E"),
        c("E", "F"), c("A", "G"), c("G", "H"), c("H", "I")
    )
}
pg_nodes <- c("R", "A", "B", "C", "D", "E", "F", "G", "H", "I")
# Two cells per node
pg_cells <- function() {
    stats::setNames(rep(pg_nodes, each = 2), sprintf("c%02d", 1:20))
}
# The expected weights of the cells of the given nodes on the path to each
# end node
pg_expected <- function(paths, cell_node = pg_cells()) {
    w <- vapply(
        paths, function(p) as.numeric(cell_node %in% p),
        numeric(length(cell_node))
    )
    rownames(w) <- names(cell_node)
    w
}

test_that("paths_from_graph() gives the cells on each root-to-end path", {
    skip_if_not_installed("igraph")
    w <- paths_from_graph(pg_cells(), pg_edges(), root = "R")
    expect_identical(w, pg_expected(list(
        D = c("R", "A", "B", "C", "D"),
        F = c("R", "A", "B", "E", "F"),
        I = c("R", "A", "G", "H", "I")
    )))
    # With the graph distance from R as pseudotime, assign_paths() labels R
    # and A as the root and B as the link of D and F
    dist <- c(
        R = 0, A = 1, B = 2, C = 3, D = 4, E = 3, F = 4, G = 2, H = 3, I = 4
    )
    cells <- assign_paths(unname(dist[pg_cells()]), w)
    expect_identical(cells$Path, rep(c(
        "root", "root", "D_links_F", "D", "D", "F", "F", "I", "I", "I"
    ), each = 2))
    expect_identical(cells$Pseudotime, unname(dist[pg_cells()]))
})

test_that("endpoints and root choose the lineages", {
    skip_if_not_installed("igraph")
    w <- paths_from_graph(pg_cells(), pg_edges(),
        root = "R", endpoints = c("I", "D")
    )
    expect_identical(w, pg_expected(list(
        I = c("R", "A", "G", "H", "I"),
        D = c("R", "A", "B", "C", "D")
    )))
    # An inner root: R becomes an end node
    w <- paths_from_graph(pg_cells(), pg_edges(), root = "B")
    expect_identical(w, pg_expected(list(
        R = c("B", "A", "R"),
        D = c("B", "C", "D"),
        F = c("B", "E", "F"),
        I = c("B", "A", "G", "H", "I")
    )))
})

test_that("edges and nodes may be given in other forms", {
    skip_if_not_installed("igraph")
    w <- paths_from_graph(pg_cells(), pg_edges(), root = "R")
    # A data frame with factor columns and an extra column of edge lengths
    edges <- data.frame(
        from = factor(pg_edges()[, 1]), to = factor(pg_edges()[, 2]),
        length = seq_len(9)
    )
    expect_identical(paths_from_graph(pg_cells(), edges, root = "R"), w)
    # The cells in another order, and nodes as a factor
    cells <- pg_cells()[20:1]
    expect_identical(
        paths_from_graph(factor(cells), pg_edges(), root = "R"),
        w[20:1, ]
    )
    # Integer node ids, compared as text
    ids <- matrix(match(pg_edges(), pg_nodes), ncol = 2)
    cells <- stats::setNames(match(pg_cells(), pg_nodes), names(pg_cells()))
    w_ids <- paths_from_graph(cells, ids, root = "1")
    expect_identical(unname(w_ids), unname(w))
    expect_identical(colnames(w_ids), c("5", "7", "10"))
})

test_that("a graph with a cycle takes the shortest path", {
    skip_if_not_installed("igraph")
    # R-A-B and R-C-D-B both reach B; the path to E goes through A
    edges <- rbind(
        c("R", "A"), c("A", "B"), c("R", "C"), c("C", "D"), c("D", "B"),
        c("B", "E")
    )
    cells <- c(x1 = "R", x2 = "A", x3 = "C", x4 = "D", x5 = "B", x6 = "E")
    w <- paths_from_graph(cells, edges, root = "R")
    expect_identical(colnames(w), "E")
    expect_identical(unname(w[, "E"]), c(1, 1, 0, 0, 1, 1))
})

test_that("paths_from_graph() checks its arguments", {
    skip_if_not_installed("igraph")
    cells <- pg_cells()
    edges <- pg_edges()
    msg_edges <- paste(
        "'edges' must be a matrix or a data frame with the two nodes of",
        "each edge in its first two columns."
    )
    msg_cells <- paste(
        "'cell_node' must be a vector of nodes named by cell, each name",
        "once."
    )
    msg_ends <- paste(
        "'endpoints' must be nodes of 'edges' other than the root,",
        "each given once."
    )
    expect_error(paths_from_graph(cells, c("R", "A"), root = "R"),
        msg_edges,
        fixed = TRUE
    )
    expect_error(paths_from_graph(cells, edges[, 1, drop = FALSE], "R"),
        msg_edges,
        fixed = TRUE
    )
    expect_error(paths_from_graph(cells, edges[0, ], "R"), msg_edges,
        fixed = TRUE
    )
    bad <- edges
    bad[2, 2] <- NA
    expect_error(paths_from_graph(cells, bad, "R"),
        "The nodes in 'edges' must not be missing or empty.",
        fixed = TRUE
    )
    expect_error(paths_from_graph(cells, edges, root = c("R", "A")),
        "'root' must be a single character string.",
        fixed = TRUE
    )
    expect_error(paths_from_graph(cells, edges, root = "Z"),
        "The root node 'Z' is not a node of 'edges'.",
        fixed = TRUE
    )
    expect_error(paths_from_graph(unname(cells), edges, "R"), msg_cells,
        fixed = TRUE
    )
    bad <- cells
    names(bad)[2] <- "c01"
    expect_error(paths_from_graph(bad, edges, "R"), msg_cells, fixed = TRUE)
    expect_error(paths_from_graph(as.list(cells), edges, "R"), msg_cells,
        fixed = TRUE
    )
    bad <- cells
    bad[3] <- "Z"
    expect_error(paths_from_graph(bad, edges, "R"),
        "Every node in 'cell_node' must be a node of 'edges'.",
        fixed = TRUE
    )
    bad[3] <- NA
    expect_error(paths_from_graph(bad, edges, "R"),
        "Every node in 'cell_node' must be a node of 'edges'.",
        fixed = TRUE
    )
    for (ends in list("R", c("D", "D"), "Z", character(0), NA_character_)) {
        expect_error(paths_from_graph(cells, edges, "R", endpoints = ends),
            msg_ends,
            fixed = TRUE
        )
    }
    # A cycle has no node with one edge
    ring <- rbind(c("R", "A"), c("A", "B"), c("B", "R"))
    expect_error(paths_from_graph(c(x = "A"), ring, "R"),
        "The graph has no end node other than the root.",
        fixed = TRUE
    )
    # A second component with end nodes X and Y
    expect_error(paths_from_graph(cells, rbind(edges, c("X", "Y")), "R"),
        "No path from the root node 'R' to the end node(s) 'X', 'Y'.",
        fixed = TRUE
    )
})
