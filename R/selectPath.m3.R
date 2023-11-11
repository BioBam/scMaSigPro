#' Subset a CDS object interactively with Shiny
#'
#' @param cdsObj An cdsObject of class `scMaSigProClass`. This cdsObject will be checked
#'   to ensure it's the right type.
#' @param redDim Dimension to use for the plot
#' @param annotation_col A character vector indicating the paths to be selected.
#' @param pseudotime_col Name of the column with Pseudotime
#' @param path_col Name of the column with Path
#'
#' @return A `scMaSigProClass` object, subsetted based on the specified paths.
#'
#' @importFrom igraph get.data.frame
#'
#' @export
#'
#'
selectPath.m3 <- function(cdsObj, redDim = "umap",
                          annotation_col = "cell.type",
                          pseudotime_col = "Pseudotime",
                          path_col = "Path") {
  # Validate is supplied opject is a valid
  assert_that(is(cdsObj, "cell_data_set"),
    msg = "Please supply a valid monocle3 cdsObject"
  )
  # Check whether the lower dimensions are calculated
  assert_that(nrow(as.data.frame(reducedDims(cdsObj)[[toupper(redDim)]])) > 1,
    msg = paste(redDim, "not found, in the cdsObj")
  )

  # Check whether the lower dimensions are calculated
  assert_that(nrow(as.data.frame(reducedDims(cdsObj)[[toupper(redDim)]])) == ncol(cdsObj),
    msg = paste("Dimensions of", redDim, "do not correspond to dimensions of counts")
  )

  # Extract UMAP
  dims <- reducedDims(cdsObj)[[toupper(redDim)]] %>% as.data.frame()

  # Set rownames as cells
  dims[["cell"]] <- rownames(dims)

  # Check whether the lower dimensions are calculated
  assert_that(all(rownames(dims) == colnames(cdsObj)),
    msg = paste("Cell Barcodes do not among", redDim, "and counts")
  )

  # Check if supplied annotation_col exist in the cdsObj
  assert_that(annotation_col %in% names(cdsObj@colData),
    msg = paste(annotation_col, "does not exist in the cell.level metadata")
  )

  # Extract the vertex cell relationships
  assert_that(!is.null(cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$pr_graph_cell_proj_closest_vertex),
    msg = paste("Vertex information is missing")
  )

  # Check pseudotime
  assert_that(!is.null(cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$pseudotime),
    msg = paste("No Pseudotime information found")
  )

  # Create Pseudotime frame
  pTime.frame <- data.frame(
    Pseudotime = cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$pseudotime,
    cell = names(cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$pseudotime)
  )

  # Create close vertex frames
  vertex.relation.frame <- data.frame(
    node = paste("Y", cdsObj@principal_graph_aux[[toupper(redDim)]]$pr_graph_cell_proj_closest_vertex[, 1], sep = "_"),
    cell = names(cdsObj@principal_graph_aux[[toupper(redDim)]]$pr_graph_cell_proj_closest_vertex[, 1])
  )

  # Create annotation_col df
  anno.df <- data.frame(
    cell = rownames(cdsObj@colData),
    anno = cdsObj@colData[[annotation_col]]
  )

  # Check before merge
  assert_that(all(anno.df[["cell"]] == dims[["cell"]]),
    msg = paste("Cells in lower dimensions does not match with cells for which annotation_col is supplied")
  )

  # Merge Anno.df with pseudotime
  anno.df <- merge(anno.df, pTime.frame, by = "cell")

  # Merge Anno.df with LD
  anno.df <- merge(anno.df, dims, by = "cell")

  # Merge with the close vertex reference
  anno.df <- merge(vertex.relation.frame, anno.df, by = "cell")

  # Set Columns
  colnames(anno.df) <- c("cell", "node", "anno", pseudotime_col, "x", "y")

  # Remove frame
  pTime.frame <- dims <- NULL

  # Extract the graph and MST
  assert_that(!is.null(cdsObj@principal_graph@listData[[toupper(redDim)]]),
    msg = paste("Principal Graph not found in cdsObj")
  )
  assert_that(!is.null(cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$dp_mst),
    msg = paste("MST not found in the cdsObj")
  )

  # Extract
  pgraph.coords <- as.data.frame(t(cdsObj@principal_graph_aux@listData[[toupper(redDim)]][["dp_mst"]]))
  pgraph.coords[["node"]] <- rownames(pgraph.coords)
  names(pgraph.coords) <- c("x", "y", "node")
  pgraph <- cdsObj@principal_graph@listData[[toupper(redDim)]]

  # Check for inf time
  if (any(is.infinite(anno.df[[pseudotime_col]]))) {
    # Check number of partitions
    numPartitions <- c(1:(length(unique(as.vector(cdsObj@clusters@listData[[toupper(redDim)]]$partitions))) - 1))

    # Create weight columns
    weight_colnames <- c(paste("weight", numPartitions, sep = "_"), "weight")
  } else {
    weight_colnames <- "weight"
  }

  # Get edge Data
  edges_df <- get.data.frame(pgraph, what = "edges")

  # Merge with Edges
  edges_df <- merge(edges_df, pgraph.coords, by.x = "from", by.y = "node")
  colnames(edges_df) <- c("from", "to", weight_colnames, "x_from", "y_from")

  # Create trajectory DF
  trajectory.df <- merge(edges_df, pgraph.coords, by.x = "to", by.y = "node")
  colnames(trajectory.df) <- c("from", "to", weight_colnames, "x_from", "y_from", "x_to", "y_to")

  # # Run Shiny
  # selection.list <- list(
  #   root = "Y_14",
  #   path1 = c("Y_2","Y_14","Y_29","Y_33","Y_34","Y_35","Y_40"),
  #   path2 = c("Y_6","Y_14","Y_16","Y_19","Y_27","Y_41")
  # )

  # View(trajectory.df)
  # View(anno.df)
  # View(pgraph.coords)
  #
  selection.list <- shinySelect(
    trajectory_data = trajectory.df,
    annotation_data = anno.df,
    label_coords = pgraph.coords,
    inputType = "Monocle3",
    pseudotime_colname = pseudotime_col
  )

  if (is.null(selection.list)) {
    warning("Nothing Returned")
  } else {
    # Create Cell Metadata
    cell.metadata <- cdsObj@colData %>% as.data.frame()

    # Subset the vertex relation frame
    vertex.relation.frame.sub <- vertex.relation.frame[vertex.relation.frame$node %in% unique(
      unlist(
        selection.list,
        use.names = FALSE
      )
    ), , drop = FALSE]

    # Subset the cell-meta
    cell.metadata.sub <- cell.metadata[rownames(cell.metadata) %in% vertex.relation.frame.sub$cell, , drop = FALSE]

    # Add annotation_col in the frame
    # cell.metadata.sub[[path_col]] <- NA

    # Extract nodes
    path1_nodes <- selection.list$path1
    path2_nodes <- selection.list$path2
    root_nodes <- selection.list$root

    # Remove roots
    path1_nodes <- path1_nodes[path1_nodes != root_nodes]
    path2_nodes <- path2_nodes[path2_nodes != root_nodes]

    # Generate subset vectors
    path1_cells <- vertex.relation.frame.sub[vertex.relation.frame.sub$node %in% path1_nodes, , drop = FALSE]
    path2_cells <- vertex.relation.frame.sub[vertex.relation.frame.sub$node %in% path2_nodes, , drop = FALSE]
    # root_cells <- vertex.relation.frame.sub[vertex.relation.frame.sub$node %in% root_nodes, , drop = F]

    # Add Root
    cell.metadata.sub[rownames(cell.metadata.sub) %in% path1_cells$cell, path_col] <- "Path1"
    cell.metadata.sub[rownames(cell.metadata.sub) %in% path2_cells$cell, path_col] <- "Path2"

    cell.metadata.sub <- cell.metadata.sub[!is.na(cell.metadata.sub[[path_col]]), , drop = FALSE]
    # cell.metadata.sub[rownames(cell.metadata.sub) %in% root_cells$cell, path_col] <- "Root"

    # Attach Pseudotime Info
    anno.df.sub <- anno.df[anno.df$cell %in% rownames(cell.metadata.sub), , drop = FALSE]
    cell.metadata.sub <- cbind(cell.metadata.sub, anno.df.sub)

    # Extract Counts
    rawCounts <- cdsObj@assays@data@listData$counts
    rawCounts <- rawCounts[, colnames(rawCounts) %in% rownames(cell.metadata.sub), drop = FALSE]

    # return(list(rawCounts=rawCounts,
    #             cell.metadata.sub=cell.metadata.sub)
    #        )
    # Call the ScMaSigPro Creator
    scmpObj <- create_scmpObj(
      counts = rawCounts,
      cell_data = cell.metadata.sub,
      pseudotime_colname = pseudotime_col,
      path_colname = path_col,
      use_as_bin = FALSE
    )
    return(scmpObj)
  }
}
