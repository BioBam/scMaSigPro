#' @title Select branching paths from a 'Cell Dataset' object from Monocle3
#'
#' @description
#' `m3_select_path()` helps select branching paths from a from a 'Cell Dataset'
#' object from Monocle3 TI Analysis. This function also has an in-built shiny app
#' thats enabled interactive selection.
#'
#' @importFrom igraph get.data.frame
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom dplyr pull arrange n
#'
#' @param cds An object from Monocle3.
#' @param latent_space Latent dimension to use for the plotting.
#' (Default is "umap")
#' @param anno_col A character string representing the column name for cell level
#' metadata containing cell level annotations. (Default is "cell_type").
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values. (Default is "Pseudotime")
#' @param path_col A character string representing the column name for branching
#' path assignment. (Default is `path_prefix`)
#' @param use_shiny Enable the selection for shiny-selection wizard.
#' (Default is TRUE)
#' @param m3_pp A list containing the character vectors for principal
#' points ("Y_") to be used for the branch selection. See examples.
#' @param plot_purity Plot a bar plot, showing the the principal points against
#' the of `anno_col`. This will show the count of annotations for each of
#' the principal points. (Default is TRUE)
#'
#' @return An object of class \code{\link{ScMaSigPro}},
#' subsetted based on the specified paths.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @export
m3_select_path <- function(cds,
                           latent_space = "umap",
                           anno_col = "cell_type",
                           ptime_col = "Pseudotime",
                           path_col = "Path",
                           use_shiny = TRUE,
                           plot_purity = FALSE,
                           m3_pp = list(
                             root_pp = c(""), path1_pp = c(""), path2_pp = c(""),
                             path1_name = NULL, path2_name = NULL
                           )) {
  if (use_shiny == FALSE) {
    # Check whether the lower dimensions are calculated
    assertthat::assert_that(
      all(names(m3_pp) %in% c(
        "root_pp", "path1_pp", "path2_pp", "path1_name",
        "path2_name"
      )),
      msg = paste("Accepted values should be one of 'root_pp',
                              'path1_pp','path2_pp', 'path1_name','path2_name'")
    )
    if (is.null(m3_pp[["path1_name"]])) {
      m3_pp[["path1_name"]] <- "Path1"
    }
    if (is.null(m3_pp[["path2_name"]])) {
      m3_pp[["path2_name"]] <- "Path2"
    }
  }

  # Set global variables
  node <- "node"
  Pseudotime <- "Pseudotime"
  median_pseudotime <- "median_pseudotime"
  anno <- "anno"
  count <- "count"

  # Validate is supplied opject is a valid
  assertthat::assert_that(is(cds, "cell_data_set"),
    msg = "Please supply a valid monocle3 cdsect"
  )
  # Check whether the lower dimensions are calculated
  assertthat::assert_that(nrow(as.data.frame(reducedDims(cds)[[toupper(latent_space)]])) > 1,
    msg = paste(latent_space, "not found, in the cds")
  )

  # Check whether the lower dimensions are calculated
  assertthat::assert_that(nrow(as.data.frame(reducedDims(cds)[[toupper(latent_space)]])) == ncol(cds),
    msg = paste("Dimensions of", latent_space, "do not correspond to dimensions of counts")
  )


  # Extract UMAP
  dims <- reducedDims(cds)[[toupper(latent_space)]] %>% as.data.frame()

  # Set rownames as cells
  dims[["cell"]] <- rownames(dims)

  # Check whether the lower dimensions are calculated
  assertthat::assert_that(all(rownames(dims) == colnames(cds)),
    msg = paste("Cell Barcodes do not among", latent_space, "and counts")
  )

  # Check if supplied anno_col exist in the cds
  assertthat::assert_that(anno_col %in% names(cds@colData),
    msg = paste(anno_col, "does not exist in the cell.level metadata")
  )

  # Extract the vertex cell relationships
  assertthat::assert_that(!is.null(cds@principal_graph_aux@listData[[toupper(latent_space)]]$pr_graph_cell_proj_closest_vertex),
    msg = paste("Vertex information is missing")
  )

  # Check pseudotime
  assertthat::assert_that(!is.null(cds@principal_graph_aux@listData[[toupper(latent_space)]][["pseudotime"]]),
    msg = paste("No Pseudotime information found")
  )

  # Create Pseudotime frame
  pTime.frame <- data.frame(
    Pseudotime = cds@principal_graph_aux@listData[[toupper(latent_space)]][["pseudotime"]],
    cell = names(cds@principal_graph_aux@listData[[toupper(latent_space)]][["pseudotime"]])
  )

  # Create close vertex frames
  vertex.relation.frame <- data.frame(
    node = paste("Y", cds@principal_graph_aux[[toupper(latent_space)]]$pr_graph_cell_proj_closest_vertex[, 1], sep = "_"),
    cell = names(cds@principal_graph_aux[[toupper(latent_space)]]$pr_graph_cell_proj_closest_vertex[, 1])
  )

  # Create anno_col df
  anno.df <- data.frame(
    cell = rownames(cds@colData),
    anno = cds@colData[[anno_col]]
  )

  # Check before merge
  assertthat::assert_that(all(anno.df[["cell"]] == dims[["cell"]]),
    msg = paste("Cells in lower dimensions does not match with cells for which anno_col is supplied")
  )

  # Merge Anno.df with pseudotime
  anno.df <- merge(anno.df, pTime.frame, by = "cell")

  # Merge Anno.df with LD
  anno.df <- merge(anno.df, dims, by = "cell")

  # Merge with the close vertex reference
  anno.df <- merge(vertex.relation.frame, anno.df, by = "cell")

  # Set Columns
  colnames(anno.df) <- c("cell", "node", "anno", ptime_col, "x", "y")

  # Remove frame
  pTime.frame <- dims <- NULL

  # Extract the graph and MST
  assertthat::assert_that(!is.null(cds@principal_graph@listData[[toupper(latent_space)]]),
    msg = paste("Principal Graph not found in cds")
  )
  assertthat::assert_that(!is.null(cds@principal_graph_aux@listData[[toupper(latent_space)]]$dp_mst),
    msg = paste("MST not found in the cds")
  )

  # Extract
  pgraph.coords <- as.data.frame(t(cds@principal_graph_aux@listData[[toupper(latent_space)]][["dp_mst"]]))
  pgraph.coords[["node"]] <- rownames(pgraph.coords)
  names(pgraph.coords) <- c("x", "y", "node")
  pgraph <- cds@principal_graph@listData[[toupper(latent_space)]]

  # Check for inf time
  if (any(is.infinite(anno.df[[ptime_col]]))) {
    # Check number of partitions
    numPartitions <- c(1:(length(unique(as.vector(cds@clusters@listData[[toupper(latent_space)]]$partitions))) - 1))

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

  # Get supplied nodes
  supplied_nodes <- unlist(m3_pp[c("root_pp", "path1_pp", "path2_pp")])

  # Run Shiny Select
  if (use_shiny) {
    selection.list <- shinySelect(
      trajectory_data = trajectory.df,
      annotation_data = anno.df,
      label_coords = pgraph.coords,
      inputType = "Monocle3",
      ptime_col = ptime_col
    )
  } else if (plot_purity && !all(supplied_nodes %in% anno.df[["node"]])) {
    # Tranfer data
    data <- anno.df

    # Order nodes by their median Pseudotime
    node_order <- data %>%
      group_by(node) %>%
      summarise(median_pseudotime = median(!!sym(Pseudotime), na.rm = TRUE)) %>%
      arrange(!!sym(median_pseudotime)) %>%
      pull(!!sym(node))

    # Convert the 'node' column to an ordered factor based on median Pseudotime
    data$node <- factor(data$node, levels = node_order)

    # Now, count the number of instances for each 'anno' within each 'node' and calculate fractions
    data_summary <- data %>%
      group_by(!!sym(node), !!sym(anno)) %>%
      summarise(count = n(), .groups = "drop")

    data_summary <- data_summary[data_summary[[count]] >= (mean(data_summary$count) + sd(data_summary[[count]])), ]

    # Plotting the data
    fraction_bar <- ggplot(data_summary, aes(x = .data$node, y = .data$count, fill = .data$anno)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(x = "Prinicpal Node", y = "Count of cells", fill = "Annotation") +
      ggtitle("Number of cells per principal point",
        subtitle = "Decide nodes for  'root_pp', 'path1_pp','path2_pp'"
      )

    return(fraction_bar)
  } else if (all(supplied_nodes %in% anno.df[["node"]])) {
    selection.list <- list(
      path1 = m3_pp[["path1_pp"]],
      path2 = m3_pp[["path2_pp"]],
      root = m3_pp[["root_pp"]],
      path1_label = m3_pp[["path1_name"]],
      path2_label = m3_pp[["path2_name"]]
    )
  }

  if (is.null(selection.list)) {
    warning("Nothing Returned")
  } else {
    # Create Cell Metadata
    cell.metadata <- cds@colData %>% as.data.frame()

    # Subset the vertex relation frame
    vertex.relation.frame.sub <- vertex.relation.frame[vertex.relation.frame$node %in% unique(
      unlist(
        selection.list,
        use.names = FALSE
      )
    ), , drop = FALSE]

    # Subset the cell-meta
    cell.metadata.sub <- cell.metadata[rownames(cell.metadata) %in% vertex.relation.frame.sub$cell, , drop = FALSE]

    # Add anno_col in the frame
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
    cell.metadata.sub[rownames(cell.metadata.sub) %in% path1_cells$cell, path_col] <- selection.list[["path1_label"]]
    cell.metadata.sub[rownames(cell.metadata.sub) %in% path2_cells$cell, path_col] <- selection.list[["path2_label"]]
    # cell.metadata.sub[rownames(cell.metadata.sub) %in% root_cells$cell, path_col] <- "Root"

    cell.metadata.sub <- cell.metadata.sub[!is.na(cell.metadata.sub[[path_col]]), , drop = FALSE]

    # Attach Pseudotime Info
    anno.df.sub <- anno.df[anno.df$cell %in% rownames(cell.metadata.sub), , drop = FALSE]
    cell.metadata.sub <- cbind(cell.metadata.sub, anno.df.sub)

    # Extract Counts
    rawCounts <- cds@assays@data@listData$counts
    rawCounts <- rawCounts[, colnames(rawCounts) %in% rownames(cell.metadata.sub), drop = FALSE]


    # Call the ScMaSigPro Creator
    scmpObj <- create_scmp(
      counts = rawCounts,
      cell_data = cell.metadata.sub,
      ptime_col = ptime_col,
      path_col = path_col,
      use_as_bin = FALSE
    )

    # Add annotation column
    scmpObj@Parameters@anno_col <- anno_col

    return(scmpObj)
  }
}
