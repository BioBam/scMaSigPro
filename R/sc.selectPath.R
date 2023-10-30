#' Select Paths for scMaSigPro
#'
#' This function selects paths to be used in scMaSigPro. It extracts
#' a sub-cdsObject based on the specified paths.
#'
#' @param cdsObj An cdsObject of class `scMaSigProClass`. This cdsObject will be checked
#'   to ensure it's the right type.
#' @param annotation A character vector indicating the paths to be selected.
#'
#' @return A `scMaSigProClass` object, subsetted based on the specified paths.
#'
#' @importFrom shiny observeEvent renderUI renderText h3 h4
#' @importFrom plotly ggplotly layout
#'  
#' @export
#'
#
selectPath.m3 <- function(cdsObj, redDim = "umap",
                       annotation = "cell.type") {
  # Validate is supplied opject is a valid
  assert_that(class(cdsObj)[1] == "cell_data_set",
    msg = "Please supply a valid monocle3 cdsObject"
  )

  # Extract UMAP
  dims <- reducedDims(cdsObj)[[toupper(redDim)]] %>% as.data.frame()

  # Extract the graph
  pgraph <- cdsObj@principal_graph@listData[[toupper(redDim)]]

  # Get edge Data
  edges_df <- get.data.frame(pgraph, what = "edges")

  # Extract Cell
  anno.df <- data.frame(
    cell = rownames(cdsObj@colData),
    anno = cdsObj@colData[[annotation]]
  )

  # Attach dim
  anno.df <- cbind(dims, anno.df)
  colnames(anno.df) <- c("x", "y", "cell", "anno")

  # Extract UMAP cordinates
  pgraph.coords <- cdsObj@principal_graph_aux@listData[[toupper(redDim)]]$dp_mst

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
      titlePanel("scMaSigPro-Monocle3"),
      sidebarLayout(
        sidebarPanel(
          h3("Interactively Select Branching Paths"),
          h5("Procedure: First select the points on in the plot and then use the following buttons to set labels. Please make sure to follow steps as directed."),
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
          actionButton("close", "Close Session"),
          actionButton("redo", "Redo Selection"),
          HTML("<hr>"),
          h4("scMaSigPro Logs"),
          verbatimTextOutput("logtext"),
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
    ),

    # Server
    server = function(input, output, session) {
      paths_selected <- reactiveVal(list())
      root_node <- reactiveVal(NULL)
      path1 <- reactiveVal(NULL)
      path2 <- reactiveVal(NULL)
      showSubTrPlot <- reactiveVal(FALSE)
      logs <- reactiveVal("")

      # Set colors
      # get_palette <- reactive({
      #     if (input$trPlotCellColor == "magma") {
      #         return("magma")
      #     } else if (input$trPlotCellColor == "nativeR") {
      #         return("nativeR")
      #     } else if (input$trPlotCellColor == "Dark2") {
      #         return("Dark2")
      #     } else if (input$trPlotCellColor == "okabe") {
      #         return("okabe")
      #     } else if (input$trPlotCellColor == "Set1") {
      #         return("set1")
      #     }
      # })

      # Main Monocle3 Plot Rendering
      output$trajectoryPlot <- renderPlotly({
        # colOption <- get_palette()

        trajectory.map <- ggplot() +
          geom_point(
            data = anno.df, aes(x = x, y = y, color = anno), size = input$trPlotCellSize, alpha = input$trPlotCellAlpha,
            stroke = input$trPlotCellStroke
          ) +
          geom_segment(data = edges_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), size = input$trPlotSegSize) +
          geom_point(data = coords_df, aes(x = x, y = y), size = input$trPlotNodeSize, color = "black") +
          geom_text(data = coords_df, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5, size = input$trPlotNodeText) +
          theme_minimal() +
          xlab("UMAP-1") +
          ylab("UMAP-2")
        #
        # if(colOption == "nativeR"){
        #     trajectory.map <- trajectory.map
        # }else if(colOption == "magma"){
        #     trajectory.map <- trajectory.map + scale_color_viridis(discrete = TRUE)
        # }

        # return
        ggplotly(trajectory.map, source = "trajectoryPlot") %>% plotly::layout(dragmode = "lasso")
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
          logs(paste(logs(), "\nSelected Root Node:", root_node()))
        } else {
          logs(paste(logs(), "\nPlease select only one node as root node."))
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
            logs(paste(logs(), "\nPlease Select root node first."))
          } else {
            # Check if exist in root
            if (root_node() %in% path1_vector) {
              # return
              path1(path1_vector)

              # Log
              logs(paste(logs(), "\nSelected nodes for Path1:", paste(path1(), collapse = ", ")))
            } else {
              logs(paste(logs(), "\nRoot node is not a part of selected Path1", root_node(), "should exist in Path1"))
            }
          }
        } else {
          logs(paste(logs(), "\nPlease Select more than two nodes to make a path."))
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
            logs(paste(logs(), "\nPlease Select root node first."))
          } else {
            # Check if exist in root
            if (root_node() %in% path2_vector) {
              # return
              path2(path2_vector)

              # Log
              logs(paste(logs(), "\nSelected nodes for Path2:", paste(path2(), collapse = ", ")))
            } else {
              logs(paste(logs(), "\nRoot node is not a part of selected Path2", root_node(), "should exist in Path2"))
            }
          }
        } else {
          logs(paste(logs(), "\nPlease Select more than two nodes to make a path."))
        }
      })


      # Subset the trajectory
      observeEvent(input$sub, {
        showSubTrPlot(TRUE)
      })

      # Log display
      output$logtext <- renderText({
        logs()
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

        # Assuming paths_selected has the paths
        path1 <- path1()
        path2 <- path2()

        # Subset the plot
        filtered_df <- coords_df[coords_df$node %in% unique(c(path1, path2)), ]
        filtered_edges <- edges_df[edges_df$from %in% unique(c(path1, path2)) & edges_df$to %in% unique(c(path1, path2)), ]

        # Get closest cells
        close.frame <- data.frame(
          nodes = cdsObj@principal_graph_aux@listData$UMAP$pr_graph_cell_proj_closest_vertex,
          cell = rownames(cdsObj@principal_graph_aux@listData$UMAP$pr_graph_cell_proj_closest_vertex)
        )
        # Attach Y
        close.frame$nodes <- paste("Y", close.frame$nodes, sep = "_")

        # Subset
        close.frame <- close.frame[close.frame$nodes %in% unique(c(path1, path2)), , drop = F]
        
        # Subset
        sub.anno.df <- anno.df[anno.df$cell %in% close.frame$cell, ]

        # Plot subset
        sub.trajectory.map <- ggplot() +
          geom_point(
            data = sub.anno.df, aes(x = x, y = y, color = anno), size = input$trPlotCellSize, alpha = input$trPlotCellAlpha,
            stroke = input$trPlotCellStroke
          ) +
          geom_segment(data = filtered_edges, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), size = input$trPlotSegSize) +
          geom_point(data = filtered_df, aes(x = x, y = y), size = input$trPlotNodeSize,color = "black") +
          geom_text(data = filtered_df, aes(x = x, y = y, label = node), vjust = 1.5, hjust = 0.5, size = input$trPlotNodeText) +
          theme_minimal() +
          xlab("UMAP-1") +
          ylab("UMAP-2")

        ggplotly(sub.trajectory.map)
      })

      # Save button behavior
      observeEvent(input$save, {
        # Get variables
        if (length(path1()) > 2 & length(path2()) > 2 & length(root_node()) == 1) {
          subset_variables <- list(root = root_node(), path1 = path1(), path2 = path2())
          print(subset_variables)
        } else {
          logs("Please Select nodes")
        }
      })
      
      observeEvent(input$close, {
          stopApp()
      })

      # Redo Selection
      observeEvent(input$redo, {
        paths_selected <- reactiveVal(list())
        root_node <- reactiveVal(NULL)
        path1 <- reactiveVal(NULL)
        path2 <- reactiveVal(NULL)
        showSubTrPlot <- reactiveVal(FALSE)
        logs <- reactiveVal("Please Select nodes")
      })
    }
  )
}
