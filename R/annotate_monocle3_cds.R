#' @title Annotate Monocle3 Object
#'
#' @description
#' `annotate_monocle3_cds()` annotates the Monocle3's cds/CellDataSet (CDS) object. It uses
#' the principal points, root cells, and minimum spanning tree to update the CDS
#' with end-point information that is used as path information.
#'
#' @param cds A Monocle3 cds/CellDataSet object (CDS).
#' @param reduction_method The dimensionality reduction method used. Default is
#' "umap". Monocle3 currently supports "UMAP" only for most procedure.
#' @param path_prefix Prefix used to annotate the paths. (Default is "Path").
#' @param root_label Label used to annotate root cells. (Default is "root").
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path. It is generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is `path_prefix`).
#' @param pseudotime_colname Name of the column in `cell.metadata` storing
#' information for Pseudotime. It is generated usingn`colData` from the
#' \pkg{SingleCellExperiment} package. (Default is "Pseudotime").
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @details
#' The function first extracts the principal points and root cells from
#' the CDS. It then calculates the minimum spanning tree and the network matrix,
#' which contains the cells along each path from root to endpoint. Finally, it
#' updates the CDS with the endpoint data and the updated cell metadata.
#'
#' @return Annotated CDS object with end-point/path information:
#'
#' @seealso
#' \code{\link[monocle3]{principal_graph_aux}}, \code{\link[monocle3]{principal_graph}},
#' \code{\link[igraph]{shortest_paths}}, `colData` from the \pkg{SingleCellExperiment} package
#'
#' @examples
#' \dontrun{
#' # Assuming 'cds' is a Monocle3 CellDataSet (CDS)
#' result <- annotate_monocle3_cds(cds,
#'   reduction_method = "umap",
#'   path_prefix = "Path",
#'   root_label = "root",
#'   path_colname = path_prefix,
#'   verbose = TRUE
#' )
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom monocle3 principal_graph principal_graph_aux pseudotime
#' @importFrom igraph degree shortest_paths
#' @importFrom SingleCellExperiment colData
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom S4Vectors DataFrame
#' @import SingleCellExperiment
#' @export
annotate_monocle3_cds <- function(cds, reduction_method = "umap",
                                  path_prefix = "Path",
                                  root_label = "root",
                                  pseudotime_colname = "Pseudotime",
                                  path_colname = path_prefix,
                                  verbose = TRUE) {
  # Convert the reduction method to upper case
  reduction_method <- toupper(reduction_method)

  # Extract principal points and convert to a data frame
  y_to_cells.df <- principal_graph_aux(cds)[[reduction_method]][["pr_graph_cell_proj_closest_vertex"]] %>%
    as.data.frame()

  # Add a new column 'barcode' with the row names and a new column 'Y' with V1 values prefixed with 'Y_'
  y_to_cells <- y_to_cells.df %>%
    mutate(barcode = rownames(y_to_cells.df), Y = paste0("Y_", y_to_cells.df$V1)) %>%
    select(-V1)

  # Extract the root cells
  root <- cds@principal_graph_aux[[reduction_method]][["root_pr_nodes"]]

  # Print the identified root points, if verbose is TRUE
  if (verbose) {
    message(paste("Identified ", root_label, " points are", paste(root, collapse = ", ")))
  }

  # Extract the minimum spanning tree
  mst <- principal_graph(cds)[[reduction_method]]

  # Get the endpoint nodes (degree 1 nodes) of the tree and remove the root nodes
  endpoints <- names(which(degree(mst) == 1))
  endpoints <- endpoints[!endpoints %in% root]

  # Verbose
  if (verbose) {
    message(paste("Identified endpoints points are", paste(endpoints, collapse = ", ")))
    message(paste("Number of potential", path_prefix, "identified are", paste(length(endpoints))))
  }

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
    message(paste("Longest", path_prefix, "is", maxPath[2], "with", maxPath[1], "cells"))
    message(paste("Shortest", path_prefix, "is", minPath[2], "with", minPath[1], "cells"))
  }

  # Update the cell dataset with the endpoint data
  endpoints <- colnames(cellWeights)

  # Extract cell metadata
  cell.meta <- as.data.frame(colData(cds))

  # Get the lineage path names
  path_names <- colnames(cellWeights)

  # Add the lineage path to each cell metadata one by one
  cell.meta$PrincipalPoints <- apply(cellWeights, 1, function(row) {
    idx <- match(1, row[match(path_names, colnames(cellWeights))])
    if (!is.na(idx)) path_names[idx] else NA
  })

  # Mark the root cells in the cell metadata
  cell.meta[rownames(cell.meta) %in% c(rownames(y_to_cells[y_to_cells$Y %in% root, ])), "PrincipalPoints"] <- root_label

  # Add pseudotime to the cell metadata
  if (all(apply(pseudotime.path.frame, 1, function(x) length(unique(x)) == 1))) {
    cell.meta[[pseudotime_colname]] <- pseudotime.path.frame[, 1]
    if (verbose) {
      message("Universal Pseudotime Detected, using monocle3")
    }
  } else {
    stop("Pseudotime is different for different paths, this is not a normal for Monocle3. Consider using 'extract_slingshot_components()'")
  }

  # Add path information
  cell.meta[[path_colname]] <- convert_to_path(cell.meta[["PrincipalPoints"]],
    path_prefix = path_prefix,
    root_label = root_label
  )

  # Update cell dataset with the updated cell metadata
  cds@colData <- DataFrame(cell.meta)

  # Return a list with updated cell dataset, y_to_cells data, endpoints, root, and number of cells per path
  return(cds)
}
