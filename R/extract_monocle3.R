#' @title Extract Monocle3 Components
#'
#' @description This function extracts the relevant components from a Monocle3 CellDataSet (CDS). These components include the
#' principal points, root cells, minimum spanning tree, network matrix, pseudotime, number of cells per path, and cell metadata.
#' It then updates the CDS with the endpoint data and the updated cell metadata.
#'
#' @param cds A Monocle3 CellDataSet (CDS).
#' @param reduction_method The dimensionality reduction method used. Default is "umap".
#' @param verbose Logical. If TRUE, print the identified root points and information about the longest and shortest paths. Default is TRUE.
#'
#' @details The function first extracts the principal points and root cells from the CDS. It then calculates the minimum spanning tree
#' and the network matrix, which contains the cells along each path from root to endpoint.
#' The function also calculates the pseudotime for each path and the number of cells per path.
#' Finally, it updates the CDS with the endpoint data and the updated cell metadata.
#'
#' @return A list containing:
#' \itemize{
#'  \item{cds}{The updated Monocle3 CellDataSet.}
#'  \item{y_to_cells}{A dataframe containing the Y coordinates of the cells.}
#'  \item{endpoints}{The endpoint nodes of the minimum spanning tree.}
#'  \item{root}{The root nodes of the minimum spanning tree.}
#'  \item{num_cells_per_path}{A dataframe containing the number of cells for each path.}
#' }
#'
#' @seealso \code{\link[monocle3]{principal_graph_aux}}, \code{\link[monocle3]{principal_graph}}, \code{\link[igraph]{shortest_paths}}, \code{\link[monocle3]{pseudotime}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'cds' is a Monocle3 CellDataSet (CDS)
#' result <- extract_monocle3_components(cds, reduction_method = "umap", verbose = TRUE)
#' }
#' 
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' 
#' @importFrom monocle3 principal_graph principal_graph_aux pseudotime
#' @importFrom igraph degree
#'
#' @export
extract_monocle3_components <- function(cds, reduction_method = "umap", verbose = TRUE) {
  # Convert the reduction method to upper case
  reduction_method <- toupper(reduction_method)

  # Extract principal points and convert to a data frame
  y_to_cells.df <- principal_graph_aux(cds)[[reduction_method]][["pr_graph_cell_proj_closest_vertex"]] %>%
    as.data.frame()
  
  # Add a new column 'barcode' with the row names and a new column 'Y' with V1 values prefixed with 'Y_'
  y_to_cells <- y_to_cells.df %>%
    mutate(barcode = rownames(y_to_cells.df), Y = paste0("Y_", y_to_cells.df$V1)) %>%
    select(-y_to_cells.df$V1)

  # Extract the root cells
  root <- cds@principal_graph_aux[[reduction_method]][["root_pr_nodes"]]

  # Print the identified root points, if verbose is TRUE
  if (verbose) {
    message(paste("Identified root points are", paste(root, collapse = ", ")))
  }

  # Extract the minimum spanning tree
  mst <- principal_graph(cds)[[reduction_method]]

  # Get the endpoint nodes (degree 1 nodes) of the tree and remove the root nodes
  endpoints <- names(which(degree(mst) == 1))
  endpoints <- endpoints[!endpoints %in% root]

  # Generate the network matrix
  cellWeights <- lapply(endpoints, function(endpoint) {
    # Find the path from root to endpoint
    path <- paste("Y", as.character(
      shortest_paths(mst, root, endpoint)$vpath[[1]]
    ),
    sep = "_"
    )

    # Get the cells along the path
    path.frame <- y_to_cells[y_to_cells$Y %in% path, ]

    # Create a binary representation for each path
    path.frame <- data.frame(weights = as.numeric(colnames(cds) %in% path.frame$barcode))
    colnames(path.frame) <- endpoint
    return(path.frame)
  }) %>%
    do.call(what = "cbind") %>%
    as.matrix()
  rownames(cellWeights) <- colnames(cds)

  # Get pseudotime for each path
  pseudotime.path.frame <- matrix(pseudotime(cds),
    ncol = ncol(cellWeights),
    nrow = ncol(cds), byrow = FALSE
  )

  # Get the number of cells for each path
  num.cells.per.path <- apply(cellWeights, 2, FUN = function(x) {
    return(length(x[x == 1]))
  }) %>% as.data.frame()
  colnames(num.cells.per.path) <- "nCells"
  num.cells.per.path$endpoint <- rownames(num.cells.per.path)

  # Print longest and shortest paths, if verbose is TRUE
  if (verbose) {
    maxPath <- num.cells.per.path[num.cells.per.path[["nCells"]] == max(num.cells.per.path[, 1]), ]
    minPath <- num.cells.per.path[num.cells.per.path[["nCells"]] == min(num.cells.per.path[, 1]), ]
    message(paste("Longest path is", maxPath[2], "with", maxPath[1], "cells"))
    message(paste("Shortest path is", minPath[2], "with", minPath[1], "cells"))
  }

  # Update the cell dataset with the endpoint data
  endpoints <- colnames(cellWeights)

  # Extract cell metadata
  cell.meta <- as.data.frame(colData(cds))

  # Get the lineage path names
  path_names <- colnames(cellWeights)

  # Add the lineage path to each cell metadata one by one
  cell.meta$Path <- apply(cellWeights, 1, function(row) {
    idx <- match(1, row[match(path_names, colnames(cellWeights))])
    if (!is.na(idx)) path_names[idx] else NA
  })

  # Mark the root cells in the cell metadata
  cell.meta[rownames(cell.meta) %in% c(rownames(y_to_cells[y_to_cells$Y %in% root, ])), "Path"] <- "root"

  # Add pseudotime to the cell metadata
  if (all(apply(pseudotime.path.frame, 1, function(x) length(unique(x)) == 1))) {
    cell.meta$Pseudotime <- pseudotime.path.frame[, 1]
    if (verbose) {
      message("Universal Pseudotime Detected, using monocel3")
    }
  } else {
    stop("Pseudotime is different for different paths, this is not a normal for monocle3. Consider using 'extract_slingshot_components()'")
  }

  # Update cell dataset with the updated cell metadata
  colData(cds) <- DataFrame(cell.meta)

  # Return a list with updated cell dataset, y_to_cells data, endpoints, root, and number of cells per path
  return(list(
    cds = cds,
    y_to_cells = y_to_cells,
    endpoints = endpoints,
    root = root,
    num_cells_per_path = num.cells.per.path
  ))
}
