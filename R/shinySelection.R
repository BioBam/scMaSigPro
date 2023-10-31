#' Function to interactively select nodes from shiny
#' 
#' @importFrom shiny observeEvent renderUI renderText h3 h4 h5 showModal fluidPage 
#' @importFrom shiny reactiveVal titlePanel sidebarPanel HTML actionButton mainPanel
#' @importFrom shiny fluidRow column uiOutput sliderInput showNotification modalDialog
#' @importFrom shiny modalButton tagList removeModal runApp shinyApp sidebarLayout stopApp
#' @importFrom plotly ggplotly layout plotlyOutput renderPlotly event_data
#' 

shinySelect <- function(trajectory_data,
                        annotation_data,
                        label_coords,
                        inputType = "Monocle3"){
    
    
    View(trajectory_data)
    View(annotation_data)
    
    # SetUI
    ui <- fluidPage(
        titlePanel("scMaSigPro"),
        sidebarLayout(
            sidebarPanel(
                h3(paste("Interactively Select Branching Paths from a", inputType, "Dataset")),
                h5("Procedure: Start by selecting points on the plot, and then use the following buttons to set labels. Please ensure that you follow the steps as instructed."),
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
                actionButton("save", "Save & Close Session"),
                actionButton("redo", "Redo Selection"),
                HTML("<hr>"),
            ),
            mainPanel(
                fluidRow(
                    column(6, h3("Monocle3-Prinicpal Graph"), plotlyOutput("trajectoryPlot")),
                    column(6, uiOutput("m3PgrapSubPlotTitle"), plotlyOutput("subTrajectoryPlot"))
                ),
                fluidRow(
                    column(3, sliderInput("trPlotNodeSize", "Node Size", min = 0.1, max = 3, value = 1)),
                    column(3, sliderInput("trPlotCellSize", "Cell Size", min = 0.1, max = 2, value = 1)),
                    column(3, sliderInput("trPlotCellAlpha", "Cell Alpha", min = 0.1, max = 1, value = 0.5)),
                    column(3, sliderInput("trPlotSegSize", "Trajectory Width", min = 0.1, max = 1, value = 0.5))
                ),
                fluidRow(
                    column(3, sliderInput("trPlotNodeText", "Node Text Size", min = 0.5, max = 3, value = 1)),
                    column(3, sliderInput("trPlotCellStroke", "Cell Stroke", min = 0.1, max = 2, value = 0.7)),
                    # column(3, selectInput("trPlotCellColor", "Cell Color",
                    #                       choices = list("Viridis Magma (Psudotime Color)" = "magma",
                    #                                      "Native R" = "nativeR",
                    #                                      "Dark2 (8 Pallets)" = "Dark2",
                    #                                      "Color Blind (Okabe Ito's)" = "okabe",
                    #                                      "Set1 (11 Pallets)" = "Set11"),
                    #                       selected = "nativeR")),
                    # column(3, sliderInput("trPlotSegSize", "Trajectory Width", min = 0.1, max = 1,value = 0.5))
                )
            )
        )
    )
    
    
    if(inputType == "Monocle3"){
        # Server
        server <- function(input, output, session) {
            selected_paths <- reactiveVal(list())
            root_node <- reactiveVal(NULL)
            path1 <- reactiveVal(NULL)
            path2 <- reactiveVal(NULL)
            showSubTrPlot <- reactiveVal(FALSE)
            
            # Main Monocle3 Plot Rendering
            output$trajectoryPlot <- renderPlotly({
                path1_highlight <- path1()
                path2_highlight <- path2()
                root_node_value <- root_node()
                
                if (!is.null(label_coords$node) && !is.null(path1_highlight)) {
                    label_coords$path1_high <- ifelse(label_coords$node %in% path1_highlight, "Yes", "No")
                } else {
                    label_coords$path1_high <- "No"
                }
                
                if (!is.null(label_coords$node) && !is.null(path2_highlight)) {
                    label_coords$path2_high <- ifelse(label_coords$node %in% path2_highlight, "Yes", "No")
                } else {
                    label_coords$path2_high <- "No"
                }
                
                if (!is.null(label_coords$node) && !is.null(root_node_value)) {
                    label_coords$root_high <- ifelse(label_coords$node == root_node_value, "Yes", "No")
                } else {
                    label_coords$root_high <- "No"
                }
                
                if (all(c("path1_high", "path2_high", "root_high") %in% colnames(label_coords))) {
                    trajectory.map <- ggplot() +
                        geom_point(
                            data = annotation_data, aes(x = x, y = y, color = anno), size = input$trPlotCellSize, alpha = input$trPlotCellAlpha,
                            stroke = input$trPlotCellStroke
                        ) +
                        geom_segment(data = trajectory_data, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), linewidth = input$trPlotSegSize) +
                        geom_point(data = label_coords, aes(x = x, y = y), size = input$trPlotNodeSize,
                                   color = ifelse(label_coords$root_high == "Yes", "red",
                                                  ifelse(label_coords$path1_high == "Yes" & label_coords$path2_high == "Yes", "purple",
                                                         ifelse(label_coords$path1_high == "Yes", "green",
                                                                ifelse(label_coords$path2_high == "Yes", "blue", "black"))))) +
                        geom_text(data = label_coords, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5, size = input$trPlotNodeText) +
                        theme_minimal() +
                        xlab("UMAP-1") +
                        ylab("UMAP-2")
                    
                    ggplotly(trajectory.map, source = "trajectoryPlot") %>% plotly::layout(dragmode = "lasso")
                }
            })
            
            
            # Set Root Node
            observeEvent(input$rootNode, {
                # Get Selected Data
                lasso_data <- event_data("plotly_selected", source = "trajectoryPlot")
                lasso_data <- lasso_data[lasso_data$curveNumber == max(lasso_data$curveNumber), , drop = F]
                
                # Get the selected point (Unique from layer-1)
                pointNumber <- as.numeric(unique(lasso_data$pointNumber)) + 1
                
                # Attach the prefix for Monocle3
                root_vector <- paste("Y", pointNumber, sep = "_")
                
                # Validate
                if (!is.null(lasso_data) && length(root_vector) == 1) {
                    # Return root vector
                    root_node(root_vector)
                    
                    # Set log message
                    showNotification(paste("Selected Root Node:", root_node()))
                } else {
                    showNotification(paste("Please select only one node as root node."))
                }
            })
            
            # Set Path-1
            observeEvent(input$Path1, {
                # Get the lasso data
                lasso_data <- event_data("plotly_selected", source = "trajectoryPlot")
                lasso_data <- lasso_data[lasso_data$curveNumber == max(lasso_data$curveNumber), , drop = F]
                
                # Get the y nodes
                pointNumber <- as.numeric(unique(lasso_data$pointNumber)) + 1
                
                # Attach Prefix
                path1_vector <- paste("Y", pointNumber, sep = "_")
                
                # Validate
                if (!is.null(lasso_data) && length(path1_vector) > 2) {
                    if (is.null(root_node())) {
                        showNotification(paste("Please Select root node first."))
                    } else {
                        # Check if exist in root
                        if (root_node() %in% path1_vector) {
                            # return
                            path1(path1_vector)
                            
                            # Log
                            showNotification(paste("Selected nodes for Path1:", paste(path1(), collapse = ", ")))
                        } else {
                            showNotification(paste("Root node is not a part of selected Path1", root_node(), "should exist in Path1"))
                        }
                    }
                } else {
                    showNotification(paste("Please Select more than two nodes to make a path."))
                }
            })
            # Set Path-2
            observeEvent(input$Path2, {
                # Get the lasso data
                lasso_data <- event_data("plotly_selected", source = "trajectoryPlot")
                lasso_data <- lasso_data[lasso_data$curveNumber == max(lasso_data$curveNumber), , drop = F]
                
                # Get the y nodes
                pointNumber <- as.numeric(unique(lasso_data$pointNumber)) + 1
                
                # Attach Prefix
                path2_vector <- paste("Y", pointNumber, sep = "_")
                
                # Validate
                if (!is.null(lasso_data) && length(path2_vector) > 2) {
                    if (is.null(root_node())) {
                        showNotification(paste("Please Select root node first."))
                    } else {
                        # Check if exist in root
                        if (root_node() %in% path2_vector) {
                            # return
                            path2(path2_vector)
                            
                            # Log
                            showNotification(paste("Selected nodes for Path2:", paste(path2(), collapse = ", ")))
                        } else {
                            showNotification(paste("Root node is not a part of selected Path2", root_node(), "should exist in Path2"))
                        }
                    }
                } else {
                    showNotification(paste("Please Select more than two nodes to make a path."))
                }
            })
            
            
            # Subset the trajectory
            observeEvent(input$sub, {
                if (length(path1()) > 2 & length(path2()) > 2 & length(root_node()) == 1) {
                    showSubTrPlot(TRUE)
                } else {
                    showNotification("You have not selected any nodes yet")
                    showSubTrPlot(FALSE)
                }
            })
            
            # Conditional Final Plot Title
            output$m3PgrapSubPlotTitle <- renderUI({
                if (showSubTrPlot()) {
                    h3("Monocle3-Principal Graph-Subset")
                }
            })
            
            # Final Plot rendering
            output$subTrajectoryPlot <- renderPlotly({
                if (!showSubTrPlot()) {
                    return(NULL)
                }
                
                path1_highlight <- path1()
                path2_highlight <- path2()
                root_node_value <- root_node()
                
                path1 <- path1()
                path2 <- path2()
                
                # Subset the plot
                trajectory_data_sub_tmp <- trajectory_data[trajectory_data$from %in% unique(c(path1, path2)), ]
                trajectory_data_sub <- trajectory_data_sub_tmp[trajectory_data_sub_tmp$to %in% unique(c(path1, path2)), ]
                
                # Subset
                annotation_data_sub <- annotation_data[annotation_data$node %in% unique(c(path1, path2)), ]
                
                # Label Coords
                label_coords_sub <- label_coords[label_coords$node %in% unique(c(path1, path2)), ]
                
                if (!is.null(label_coords_sub$node)) {
                    label_coords_sub$path1_high <- ifelse(label_coords_sub$node %in% path1_highlight, "Yes", "No")
                    label_coords_sub$path2_high <- ifelse(label_coords_sub$node %in% path2_highlight, "Yes", "No")
                    label_coords_sub$root_high <- ifelse(label_coords_sub$node == root_node_value, "Yes", "No")
                }
                
                sub.trajectory.map <- ggplot() +
                    geom_point(
                        data = annotation_data_sub, aes(x = x, y = y, color = anno), size = input$trPlotCellSize, alpha = input$trPlotCellAlpha,
                        stroke = input$trPlotCellStroke
                    ) +
                    geom_segment(data = trajectory_data_sub, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), linewidth = input$trPlotSegSize) +
                    geom_point(data = label_coords_sub, aes(x = x, y = y), size = input$trPlotNodeSize,
                               color = ifelse(label_coords_sub$root_high == "Yes", "red",
                                              ifelse(label_coords_sub$path1_high == "Yes" & label_coords_sub$path2_high == "Yes", "purple",
                                                     ifelse(label_coords_sub$path1_high == "Yes", "green",
                                                            ifelse(label_coords_sub$path2_high == "Yes", "blue", "black"))))) +
                    geom_text(data = label_coords_sub, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5, size = input$trPlotNodeText) +
                    theme_minimal() +
                    xlab("UMAP-1") +
                    ylab("UMAP-2")
                
                ggplotly(sub.trajectory.map)
            })
            
            # Save button behavior
            observeEvent(input$save, {
                # Get variables
                if (length(path1()) > 2 & length(path2()) > 2 & length(root_node()) == 1) {
                    # Set Reactive
                    selected_paths(list(root = root_node(), path1 = path1(), path2 = path2()))
                    
                    # Check before close
                    if (length(selected_paths()) != 3) {
                        showModal(
                            modalDialog(
                                title = "Warning",
                                "You haven't selected paths, do you want to close the session?",
                                footer = tagList(
                                    actionButton("confirmClose", "Yes, close it!"),
                                    modalButton("No")
                                )
                            )
                        )
                    } else {
                        stopApp(selected_paths())
                    }
                } else {
                    # Check before close
                    if (length(selected_paths()) != 3) {
                        showModal(
                            modalDialog(
                                title = "Warning",
                                "You haven't selected paths, do you want to close the session?",
                                footer = tagList(
                                    actionButton("confirmClose", "Yes, close it!"),
                                    modalButton("No")
                                )
                            )
                        )
                    } else {
                        stopApp(selected_paths())
                    }
                }
            })
            
            observeEvent(input$confirmClose, {
                removeModal()
                stopApp()
            })
            
            # Redo Selection
            observeEvent(input$redo, {
                selected_paths <- reactiveVal(list())
                root_node <- reactiveVal(NULL)
                path1 <- reactiveVal(NULL)
                path2 <- reactiveVal(NULL)
                showSubTrPlot <- reactiveVal(FALSE)
            })
        }
        
    }else{
        stop("Currently only monocle3 is supported")
    }
    
    selection <- suppressMessages(runApp(shinyApp(ui, server)))
    
    if(is.null(selection)){
        return(NULL)
    }else{
        return(selection)
    }
    
} 