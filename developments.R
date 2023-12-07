library(tidyverse)
library(assertthat)
library(scMaSigPro)
library(shiny)
library(igraph)
library(plotly)

#cds <- readRDS("/supp_data/Analysis_Public_Data/rep1/rep1_cds.RDS")



ob <- m3_select_path(cds,
              redDim = "umap",
              annotation_col = "cell_type",
              path_col = "Path",
              pseudotime_col = "Pseudotime",
              use_shiny = F)
ob



pltM3VertexPurity <- function(cdsObj, redDim = "umap",
                        annotation_col = "cell_type",
                        pseudotime_col = "Pseudotime"){
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
    
    data <- anno.df
    
    # Order nodes by their median Pseudotime
    node_order <- data %>%
        group_by(node) %>%
        summarise(median_pseudotime = median(Pseudotime, na.rm = TRUE)) %>%
        arrange(median_pseudotime) %>%
        pull(node)
    
    # Convert the 'node' column to an ordered factor based on median Pseudotime
    data$node <- factor(data$node, levels = node_order)
    
    # Now, count the number of instances for each 'anno' within each 'node' and calculate fractions
    data_summary <- data %>%
        group_by(node, anno) %>%
        summarise(count = n(), .groups = 'drop') %>%
        mutate(fraction = count)
    
    data_summary <- data_summary[data_summary$fraction >= (mean(data_summary$fraction) + sd(data_summary$fraction)),]
    
    # Plotting the data
    ggplot(data_summary, aes(x = node, y = fraction, fill = anno)) +
        geom_bar(stat = 'identity', position = 'stack') +
        theme_minimal() + coord_flip()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate X labels for better readability
        labs(x = "Node", y = "Fraction", fill = "Annotation") +
        ggtitle("Fraction of Annotations per Node Ordered by Pseudotime")
    
}


pltM3VertexPurity(cdsObj = cds)
