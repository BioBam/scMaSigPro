#' Subset a CDS object interactively with Shiny
#'
#' @param cdsObj An cdsObject of class `scMaSigProClass`. This cdsObject will be checked
#'   to ensure it's the right type.
#' @param redDim Dimension to use for the plot
#' @param annotation A character vector indicating the paths to be selected.
#' @param returnType 'scmpObj' or subseted 'cds'
#'
#' @return A `scMaSigProClass` object, subsetted based on the specified paths.
#'
#' @importFrom igraph get.data.frame 
#'
#' @export
#'
selectPath.m3 <- function(cdsObj, redDim = "umap",
                          annotation = "cell.type",
                          returnType = "scmpObj") {
    
  # Validate is supplied opject is a valid
  assert_that(class(cdsObj)[1] == "cell_data_set",
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
  
  # Check if supplied annotation exist in the cdsObj
  assert_that(annotation %in% names(cds@colData),
              msg = paste(annotation, "does not exist in the cell.level metadata")
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
  pTime.frame <- data.frame(pTime = cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$pseudotime,
                            cell = names(cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$pseudotime))
  
  # Create close vertex frames
  vertex.relation.frame <- data.frame(
      node = paste("Y", cds@principal_graph_aux[[toupper(redDim)]]$pr_graph_cell_proj_closest_vertex[,1], sep = "_"),
      cell = names(cds@principal_graph_aux[[toupper(redDim)]]$pr_graph_cell_proj_closest_vertex[,1])
  )
  
  # Create Annotation df
  anno.df <- data.frame(
      cell = rownames(cdsObj@colData),
      anno = cdsObj@colData[[annotation]]
  )
  
  # Check before merge
  assert_that(all(anno.df[["cell"]] == dims[["cell"]]),
              msg = paste("Cells in lower dimensions does not match with cells for which annotation is supplied")
  )
  
  # Merge Anno.df with pseudotime
  anno.df <-  merge(anno.df, pTime.frame, by = "cell")
  
  # Merge Anno.df with LD
  anno.df <- merge(anno.df, dims, by = "cell")
  
  # Merge with the close vertex reference
  anno.df <- merge(vertex.relation.frame, anno.df, by = "cell")
  
  # Set Columns
  colnames(anno.df) <- c("cell", "node", "anno", "pTime", "x", "y")
  
  # Remove frame
  pTime.frame <- vertex.relation.frame <- dims <- NULL
  
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
  
  # Get edge Data
  edges_df <- get.data.frame(pgraph, what = "edges")
  
  # Merge with Edges
  edges_df <- merge(edges_df, pgraph.coords, by.x = "from", by.y = "node")
  colnames(edges_df) <- c("from", "to", "weight", "x_from", "y_from")
  trajectory.df <- merge(edges_df, pgraph.coords, by.x = "to", by.y = "node")
  colnames(trajectory.df) <- c("from", "to", "weight", "x_from", "y_from", "x_to", "y_to")
  
  # RUn Shiny
  selection.list <- shinySelect(trajectory_data = trajectory.df, 
              annotation_data = anno.df,
              label_coords = pgraph.coords, 
              inputType = "Monocle3")
  
  
  if(is.null(selection.list)){
      warning("Nothing Returned")
  }else{
      return(selection.list)
  }
}