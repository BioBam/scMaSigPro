#' Select Paths for scMaSigPro
#'
#' This function selects paths to be used in scMaSigPro. It extracts
#' a sub-object based on the specified paths.
#'
#' @param obj An object of class `scMaSigProClass`. This object will be checked
#'   to ensure it's the right type.
#' @param sel.path A character vector indicating the paths to be selected.
#' @param balance_paths Logical indicating if paths should be balanced. Default is TRUE.
#' @param pathCol The column name in `colData(sceObj)` that indicates the path.
#'   Defaults to "Path".
#' @param plot_paths Logical indicating if paths should be plotted. Default is TRUE.
#' @param pTimeCol The column name indicating the pseudotime. Default is "Pseudotime".
#' @param verbose Logical indicating if verbose messages should be displayed. Default is TRUE.
#'
#' @return A `SingleCellExperiment` object subsetted based on the specified paths.
#'
#' @examples
#' # Assuming you have an example object of class `scMaSigProClass`:
#' # selected_obj <- selectPath(example_obj, sel.path = c("Path1", "Path3"))
#'
#' @export
#' 
# 

selectPath <- function(obj, redDim = "umap") {

    # Validate is supplied opject is a valid
    assert_that(class(obj)[1] == "cell_data_set",
                msg = "Please supply a valid monocle3 object")

    # Extract UMAP
    dims <- reducedDims(obj)[[toupper(redDim)]] %>% as.data.frame()

    # Extract the graph
    pgraph <- obj@principal_graph@listData[[toupper(redDim)]]

    # Get edge Data
    edges_df <- get.data.frame(pgraph, what = "edges")

    # Extract UMAP cordinates
    pgraph.coords <- obj@principal_graph_aux@listData[[toupper(redDim)]]$dp_mst

    # Get the x-y
    coords_df <- as.data.frame(t(pgraph.coords))
    coords_df$node <- rownames(coords_df)
    names(coords_df) <- c("x", "y", "node") # Rename columns for clarity

    # Join the coordinates to the edges data frame
    edges_df <- merge(edges_df, coords_df, by.x = "from", by.y = "node")
    colnames(edges_df) <- c("from", "to", "weight", "x_from", "y_from")
    edges_df <- merge(edges_df, coords_df, by.x = "to", by.y = "node")
    colnames(edges_df) <- c("from", "to", "weight", "x_from", "y_from", "x_to", "y_to")

    shinyApp(
        ui = fluidPage(
            titlePanel("scMaSigPro"),
            sidebarLayout(
                sidebarPanel(
                    h2("Select Branching Paths"),
                    h3("Interactively select the two branching paths, to examine differential expression trends with scMaSigPro."),
                    h5("Procedure: First select the points on in the plot and then use the following buttons to set labels."),
                    HTML("<hr>"), 
                    h4("Step-1: Identify a root node."),
                    actionButton("rootNode", "Set as Root Node"),
                    HTML("<br>"), 
                    h4("Step-2: Identify branching Path-1"),
                    actionButton("Path1", "Set as Path-1"),
                    HTML("<br>"), 
                    h4("Step-3: Identify branching Path-2"),
                    actionButton("Path2", "Set as Path-2"),
                    HTML("<br>"), 
                    h4("Step-4: Subset the trajectory"),
                    actionButton("sub", "Subset"),
                    h4("Step-5: Save Selection"),
                    actionButton("save", "Save"),
                    actionButton("redo", "Redo Selection"),
                    HTML("<hr>"), 
                    h4("scMaSigPro Logs"),
                    verbatimTextOutput("logtext"),
                    HTML("<hr>"), 
                    h4("Additional Information"),
                ),
                mainPanel(
                    
                    fluidRow(column(6, h3("Monocle3-Prinicpal Graph"),plotlyOutput("mainPlot"))),
                    fluidRow(column(6,uiOutput("finalTitle"),plotlyOutput("finalPlot")))
                )
            )
        )
        ,

        server = function(input, output, session) {
            
            paths_selected <- reactiveVal(list())
            root_node <- reactiveVal(NULL)
            path1 <- reactiveVal(NULL)
            path2 <- reactiveVal(NULL)
            show_final_plot <- reactiveVal(FALSE)
            logs <- reactiveVal("")
            
            # Main Monocle3 Plot Rendering
            output$mainPlot <- renderPlotly({
                p <- ggplot() +
                    geom_segment(data = edges_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), size = 0.5) +
                    geom_point(data = coords_df, aes(x = x, y = y), size = 1.5, color = "black") +
                    geom_text(data = coords_df, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5) +
                    theme_minimal() + xlab("UMAP-1") + ylab("UMAP-2")
                ggplotly(p, source = "mainPlot") %>% layout(dragmode = "lasso")
            })
            
            # Set Root Node
            observeEvent(input$rootNode, {
                lasso_data <- event_data("plotly_selected", source = "mainPlot")
                pointNumber <- unique(lasso_data$pointNumber)
                if (!is.null(lasso_data) && length(pointNumber) == 1) {
                    selected_nodes <- pointNumber + 1
                    root_node(paste(coords_df$node[selected_nodes], collapse = ", "))
                    logs(paste(logs(), "\nSelected Root Node:", root_node()))
                } else {
                    logs(paste(logs(), "\nPlease Select only one node as root node.")) 
                }
            })
            
            # Set Path-1
            observeEvent(input$Path1, {
                lasso_data <- event_data("plotly_selected", source = "mainPlot")
                pointNumber <- unique(lasso_data$pointNumber)
                if (!is.null(lasso_data) && length(pointNumber) > 2) {
                    if (is.null(root_node())){
                        logs(paste(logs(), "\nPlease Select root node first.")) 
                    }else{
                        
                        # Select nodes
                        selected_nodes <- pointNumber + 1
                        
                        if (root_node() %in% coords_df$node[selected_nodes]) {
                            path1(paste(coords_df$node[selected_nodes], collapse = ", "))
                            logs(paste(logs(), "\nSelected nodes for Path1:", path1()))
                            
                            # logs(paste(logs(), "\n", path1))
                            # 
                            # if (length(paths_selected()) == 2) {
                            #     show_final_plot(TRUE)
                            # }
                        }else{
                            logs(paste(logs(), "\nRoot node is not a part of selected Path1", root_node(), "should exist in Path1"))
                        }
                    }
                } else {
                    logs(paste(logs(), "\nPlease Select more than two nodes to make a path.")) 
                }
            })
            # Set Path-2
            
            observeEvent(input$Path2, {
                lasso_data <- event_data("plotly_selected", source = "mainPlot")
                pointNumber <- unique(lasso_data$pointNumber)
                if (!is.null(lasso_data) && length(pointNumber) > 2) {
                    if (is.null(root_node())){
                        logs(paste(logs(), "\nPlease Select root node first.")) 
                    }else{
                        if (is.null(path1())){
                            logs(paste(logs(), "\nPlease select Path1 first.")) 
                        }else{
                            
                            # Select nodes
                            selected_nodes <- pointNumber + 1
                            sel_nodes <- coords_df$node[selected_nodes]
                            
                            
                            if (root_node() %in% sel_nodes) {
                                
                                print(path1())
                                print(sel_nodes)
                                
                                if(path1() == sel_nodes){
                                    logs(paste(logs(), "\nPath2 contains same nodes as path. (", path1()[path1() %in% sel_nodes], ") are common nodes"))
                                }else{
                                path2(paste("Selected Nodes:", paste(sel_nodes, collapse = ", ")))
                                logs(paste(logs(), "\nSelected nodes for Path2:", path2()))
                                }
                            }else{
                                logs(paste(logs(), "\nRoot node is not a part of selected Path2", root_node(), "should exist in Path2"))
                            }
                        }
                    }
                } else {
                    logs(paste(logs(), "\nPlease Select more than two nodes to make a path.")) 
                }
            })
            
            
            # Subset the trajectory
            observeEvent(input$sub, {
                show_final_plot(TRUE)
            })
            
            # Log display
            output$logtext <- renderText({
                logs()
            })
            
            # Conditional Final Plot Title
            output$finalTitle <- renderUI({
                if (show_final_plot()) {
                    h3("Monocle3-Principal Graph-Subset")
                }
            })
            
            # Final Plot rendering
            output$finalPlot <- renderPlotly({
                if (!show_final_plot()) return(NULL)
                
                # Assuming paths_selected has the paths
                path1 <- strsplit(paths_selected()[[1]], "Selected Nodes: ")[[1]][2]
                path2 <- strsplit(paths_selected()[[2]], "Selected Nodes: ")[[1]][2]
                nodes <- unlist(c(strsplit(path1, ", "), strsplit(path2, ", ")))
                
                filtered_df <- coords_df[coords_df$node %in% nodes, ]
                filtered_edges <- edges_df[edges_df$from %in% nodes & edges_df$to %in% nodes, ]
                
                p <- ggplot() +
                    geom_segment(data = filtered_edges, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), size = 0.5) +
                    geom_point(data = filtered_df, aes(x = x, y = y), size = 1.5, color = "black") +
                    geom_text(data = filtered_df, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5) +
                    theme_minimal()
                
                ggplotly(p)
            })
            
            # Save button behavior
            observeEvent(input$save, {
                if (length(paths_selected()) != 2) return()
                
                saved_data <- list(root = root_node(), path1 = paths_selected()[[1]], path2 = paths_selected()[[2]])
                # Here, you can perform operations to save the data, or print the saved_data list to see the saved nodes.
                print(saved_data) 
            })
            
            # Redo Selection
            observeEvent(input$redo, {
                paths_selected(list())
                root_node(NULL)
                current_path(NULL)
                show_final_plot(FALSE)
                logs("")
            })
    
        })}
    #     server = function(input, output, session) {
    #         paths_selected <- reactiveVal(list())
    #         current_path <- reactiveVal(NULL)
    #         show_final_plot <- reactiveVal(FALSE)
    # 
    #         output$mainPlot <- renderPlotly({
    #             p <- ggplot() +
    #                 geom_segment(data = edges_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), size = 0.5) +
    #                 geom_point(data = coords_df, aes(x = x, y = y), size = 1.5, color = "black") +
    #                 geom_text(data = coords_df, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5) +
    #                 theme_minimal() + xlab("UMAP-1") + ylab("UMAP-2")
    # 
    #             ggplotly(p, source = "mainPlot") %>% layout(dragmode = "lasso")
    #         })
    # 
    #         observeEvent(input$Path1, {
    #             if (is.null(current_path())) return()
    # 
    #             paths_selected(append(paths_selected(), list(current_path())))
    #             current_path(NULL)
    # 
    #             if (length(paths_selected()) == 2) {
    #                 show_final_plot(TRUE)
    #             }
    #         })
    # 
    #         observe({
    #             lasso_data <- event_data("plotly_selected", source = "mainPlot")
    #             if (!is.null(lasso_data)) {
    #                 selected_nodes <- lasso_data$pointNumber + 1
    #                 current_path(paste("Selected Nodes:", paste(coords_df$node[selected_nodes], collapse = ", ")))
    #             }
    #         })
    # 
    #         output$logtext <- renderText({
    #             current_path()
    #         })
    #         
    #         output$finalTitle <- renderUI({
    #             if (show_final_plot()) {
    #                 h3("Monocle3-Principal Graph-Subset")
    #             }
    #         })
    #         
    # 
    #         output$finalPlot <- renderPlotly({
    #             if (!show_final_plot()) return(NULL)
    # 
    #             path1 <- strsplit(paths_selected()[[1]], "Selected Nodes: ")[[1]][2]
    #             path2 <- strsplit(paths_selected()[[2]], "Selected Nodes: ")[[1]][2]
    #             nodes <- unlist(c(strsplit(path1, ", "), strsplit(path2, ", ")))
    # 
    #             filtered_df <- coords_df[coords_df$node %in% nodes, ]
    #             filtered_edges <- edges_df[edges_df$from %in% nodes & edges_df$to %in% nodes, ]
    # 
    #             p <- ggplot() +
    #                 geom_segment(data = filtered_edges, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), size = 0.5) +
    #                 geom_point(data = filtered_df, aes(x = x, y = y), size = 1.5, color = "black") +
    #                 geom_text(data = filtered_df, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5) +
    #                 theme_minimal()
    # 
    #             ggplotly(p)
    #         })
    #     }
    # )
    
    

    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # # Extract the sce class
    # sceObj <- obj@sce
    # 
    # # Extract Cell MetaData
    # cell.meta.raw <- as.data.frame(colData(sceObj))
    # 
    # # Extract Paths from the metadata
    # cell.meta.sub <- cell.meta.raw[cell.meta.raw[[pathCol]] %in% sel.path, ]
    # 
    # # Plot
    # if (plot_paths) {
    #     before.plt <- ggplot(cell.meta.sub, aes(x = .data[[pathCol]], y = .data[[pTimeCol]])) +
    #         geom_boxplot(color = "#f58a53") + # This creates the boxplots for each category
    #         stat_summary(fun = median, geom = "line", aes(group = 1), color = "#e84258") +
    #         stat_summary(fun = median, geom = "point", color = "#159287") +
    #         labs(
    #             title = "A. Before Balance",
    #             x = "Available Paths", y = "Inferred Pseudotime"
    #         ) +
    #         theme_minimal(base_size = 15)
    # }
    # 
    # # Balance
    # if (balance_paths) {
    #     # Get the avaibale paths
    #     avail_paths <- unique(cell.meta.sub[[pathCol]])
    #     
    #     # Store original ranges
    #     original_ranges <- lapply(avail_paths, function(path) {
    #         return(
    #             range(cell.meta.sub[cell.meta.sub[[pathCol]] == path, pTimeCol], na.rm = TRUE)
    #         )
    #     })
    #     names(original_ranges) <- avail_paths
    #     
    #     # Get the span of the ranges
    #     range_spans <- sapply(original_ranges, diff)
    #     
    #     # Order the ranges
    #     sorted_paths <- names(range_spans[order(range_spans)])
    #     
    #     # First path is the smallest
    #     smallest_path <- sorted_paths[1]
    #     
    #     if (!is.na(smallest_path)) {
    #         # Get the range of the smallest path
    #         smallest_range <- original_ranges[[smallest_path]]
    #         
    #         # Set rownames to columns
    #         cell.meta.sub$row_id <- rownames(cell.meta.sub)
    #         
    #         # Get the list of cells to drop
    #         drop_cells <- rownames(cell.meta.sub[!(cell.meta.sub[[pTimeCol]] >= smallest_range[1] &
    #                                                    cell.meta.sub[[pTimeCol]] <= smallest_range[2]), ])
    #         
    #         # Subset
    #         cell.meta.sub.sliced <- cell.meta.sub[!(cell.meta.sub$row_id %in% drop_cells), , drop = F]
    #         
    #         # Set the rownames
    #         rownames(cell.meta.sub.sliced) <- cell.meta.sub.sliced$row_id
    #         cell.meta.sub.sliced$row_id <- NULL
    #         
    #         if (verbose) {
    #             removed_count <- nrow(cell.meta.sub) - nrow(cell.meta.sub.sliced)
    #             message(paste(removed_count, "cells were removed to match the pseudotime range of", smallest_path))
    #         }
    #         
    #         if (plot_paths) {
    #             after.plt <- ggplot(cell.meta.sub.sliced, aes(x = .data[[pathCol]], y = .data[[pTimeCol]])) +
    #                 geom_boxplot(color = "#e84258") + # This creates the boxplots for each category
    #                 stat_summary(fun = median, geom = "line", aes(group = 1), color = "#f58a53") +
    #                 stat_summary(fun = median, geom = "point", color = "#159287") +
    #                 labs(
    #                     title = "B. After Balance",
    #                     x = "Available Paths", y = "Inferred Pseudotime"
    #                 ) +
    #                 theme_minimal(base_size = 15)
    #         }
    #         
    #         
    #         # Select Cells
    #         sceObj_sub <- sceObj[, rownames(colData(sceObj)) %in% rownames(cell.meta.sub.sliced)]
    #         
    #         # Add the Object Back
    #         obj@sce <- sceObj_sub
    #         
    #         # Plot
    #         if (plot_paths) {
    #             comb.plt <- ggarrange(before.plt, after.plt)
    #             print(comb.plt)
    #         }
    #         
    #         # return
    #         return(obj)
    #     } else {
    #         message("Nothing to remove, paths correspond")
    #         return(obj)
    #     }
    # }
#}