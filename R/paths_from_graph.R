#' @title Lineage weights from a trajectory graph
#'
#' @description
#' `paths_from_graph()` gives the lineage weights of the cells of a
#' trajectory graph, such as the principal graph of Monocle3 or a minimum
#' spanning tree of clusters, for [assign_paths()]. Each lineage is the
#' shortest path from the root node to an end node, and a cell is on it
#' (weight 1) when its node is on that path.
#'
#' @param cell_node A vector with the node of each cell, named by cell. The
#' nodes are compared with the node names of `edges` as character strings.
#' @param edges A matrix or data frame whose first two columns hold the two
#' nodes of each edge of the graph, which is taken as undirected.
#' @param root The name of the root node.
#' @param endpoints The names of the end nodes, one lineage each; the
#' default, `NULL`, takes every node with one edge other than the root, in
#' the order in which they first appear in `edges`.
#'
#' @details
#' The paths are found with `igraph::shortest_paths()` on the unweighted
#' graph, so this function needs the suggested package \pkg{igraph}. In a
#' tree, as principal graphs and minimum spanning trees are, the path from
#' the root to each end node is the only one. Every node of `cell_node`
#' must be a node of `edges`, and every end node must be reachable from the
#' root.
#'
#' @return A numeric matrix with one row per cell, in the order and with the
#' names of `cell_node`, and one column per end node, named by it, holding 1
#' when the node of the cell is on the path from the root to that end node
#' and 0 otherwise.
#'
#' @seealso [assign_paths()], which turns these weights into path labels.
#'
#' @examples
#' # A tree: the root R, a segment A-B shared by the end nodes C and D, and
#' # the end node E that branches off at A
#' edges <- rbind(
#'     c("R", "A"), c("A", "B"), c("B", "C"), c("B", "D"), c("A", "E")
#' )
#' cell_node <- c(
#'     cell1 = "R", cell2 = "A", cell3 = "B", cell4 = "C", cell5 = "D",
#'     cell6 = "E"
#' )
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'     weights <- paths_from_graph(cell_node, edges, root = "R")
#'     print(weights)
#'     # With a pseudotime that every lineage shares
#'     assign_paths(c(0, 1, 2, 3, 3, 2), weights)
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @family input helpers
#'
#' @references
#' Van den Berge K, Roux de Bézieux H, Street K, Saelens W, Cannoodt R,
#' Saeys Y, Dudoit S, Clement L (2020). Trajectory-based differential
#' expression analysis for single-cell sequencing data. Nature
#' Communications 11, 1201. \doi{10.1038/s41467-020-14766-3}
#'
#' Csárdi G, Nepusz T (2006). The igraph software package for complex
#' network research. InterJournal Complex Systems, 1695.
#' <https://igraph.org>
#'
#' @export
paths_from_graph <- function(cell_node, edges, root, endpoints = NULL) {
    .check(requireNamespace("igraph", quietly = TRUE),
        msg = paste(
            "paths_from_graph() needs the 'igraph' package; install it with",
            "install.packages('igraph')."
        )
    )
    # The graph: one row per edge, the two nodes in the first two columns
    .check(
        (is.matrix(edges) || is.data.frame(edges)) && ncol(edges) >= 2L &&
            nrow(edges) >= 1L,
        msg = paste(
            "'edges' must be a matrix or a data frame with the two nodes of",
            "each edge in its first two columns."
        )
    )
    edge_list <- cbind(as.character(edges[, 1]), as.character(edges[, 2]))
    .check(!anyNA(edge_list) && all(nzchar(edge_list)),
        msg = "The nodes in 'edges' must not be missing or empty."
    )
    nodes <- unique(as.vector(t(edge_list)))
    .check_strings(root = root)
    .check(root %in% nodes,
        msg = paste0("The root node '", root, "' is not a node of 'edges'.")
    )

    # The cells and their nodes
    cells <- names(cell_node)
    .check(
        is.atomic(cell_node) && length(cell_node) >= 1L &&
            .unique_names(cells),
        msg = paste(
            "'cell_node' must be a vector of nodes named by cell, each name",
            "once."
        )
    )
    cell_node <- as.character(cell_node)
    .check(all(cell_node %in% nodes),
        msg = "Every node in 'cell_node' must be a node of 'edges'."
    )

    # The end nodes: by default every node with one edge, except the root
    if (is.null(endpoints)) {
        degree <- table(factor(edge_list, levels = nodes))
        endpoints <- setdiff(nodes[degree == 1L], root)
        .check(length(endpoints) >= 1L,
            msg = "The graph has no end node other than the root."
        )
    } else {
        .check(
            is.character(endpoints) && length(endpoints) >= 1L &&
                !anyNA(endpoints) && !anyDuplicated(endpoints) &&
                all(endpoints %in% nodes) && !root %in% endpoints,
            msg = paste(
                "'endpoints' must be nodes of 'edges' other than the root,",
                "each given once."
            )
        )
    }

    graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    # An end node that cannot be reached is at an infinite distance from
    # the root; such nodes are reported before the paths are searched
    reach <- igraph::distances(graph, v = root, to = endpoints)
    unreached <- endpoints[!is.finite(reach[1L, ])]
    .check(length(unreached) == 0L,
        msg = paste0(
            "No path from the root node '", root, "' to the end node(s) ",
            paste0("'", unreached, "'", collapse = ", "), "."
        )
    )
    paths <- igraph::shortest_paths(graph, from = root, to = endpoints)$vpath

    weights <- matrix(0,
        nrow = length(cells), ncol = length(endpoints),
        dimnames = list(cells, endpoints)
    )
    for (j in seq_along(endpoints)) {
        weights[, j] <- as.numeric(cell_node %in% igraph::as_ids(paths[[j]]))
    }
    weights
}
