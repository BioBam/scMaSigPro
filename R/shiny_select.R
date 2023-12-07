#' Function to interactively select nodes from shiny
#'
#' @param trajectory_data Dataframe for trajectory
#' @param annotation_data Dataframe for cell UMAP
#' @param label_coords Datatframe for labels
#' @param inputType Input type
#'
#' @importFrom shiny observeEvent renderUI renderText h3 h4 h5 showModal fluidPage
#' @importFrom shiny reactiveVal titlePanel sidebarPanel HTML actionButton mainPanel
#' @importFrom shiny fluidRow column uiOutput sliderInput showNotification modalDialog
#' @importFrom shiny modalButton tagList removeModal runApp shinyApp sidebarLayout stopApp
#' @importFrom shiny verbatimTextOutput radioButtons
#' @importFrom plotly ggplotly layout plotlyOutput renderPlotly event_data config
#'
#' @keywords internal
#'

shiny_select <- function(trajectory_data,
                        annotation_data,
                        label_coords,
                        pseudotime_colname,
                        inputType = "Monocle3") {
    
    # Set viridis
    viridis_colors <- rev(c("#f0f921", "#f89540", "#cc4778", "#7e03a8", "#0d0887"))
    cell_colors <- c(
        "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
        "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
        "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
        "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
    )
    
    
  # SetUI
  ui <- fluidPage(
    titlePanel("scMaSigPro"),
    sidebarLayout(
      sidebarPanel(
        h4(paste("Interactively Select Branching Paths from a", inputType, "Dataset")),
        # h5("Please ensure that you follow the steps as instructed."),
        h4("Step-1: Identify a root node."),
        h5("Select only one node as the root."),
        actionButton("rootNode", "Set as Root Node"),
        HTML("<br>"),
        h4("Step-2: Identify branching Path-1"),
        h5("Select atleast 3 nodes (including the root node)"),
        actionButton("Path1", "Set as Path-1"),
        HTML("<br>"),
        h4("Step-3: Identify branching Path-2"),
        h5("Select atleast 3 nodes (including the root node, excluding path-1 nodes)"),
        actionButton("Path2", "Set as Path-2"),
        HTML("<br>"),
        h4("Step-4: Subset the trajectory"),
        h5("Validate before save"),
        actionButton("sub", "Subset"),
        hr(),
        h5("Set labels"),
        textInput("path1_label", label = "Path1 Label:", value = "Path1"),
        textInput("path2_label", label = "Path2 Label:", value = "Path2")
      ),
      mainPanel(
        fluidRow(
          column(6),
          column(3, actionButton("save", "Save & Close Session")),
          column(3, actionButton("redo", "Redo Selection"))
        ),
        fluidRow(
          column(6, h3("Monocle3-Prinicpal Graph"), plotlyOutput("trajectoryPlot")),
          column(6, uiOutput("m3PgrapSubPlotTitle"), plotlyOutput("subTrajectoryPlot"))
        ),
        fluidRow(
          column(4, verbatimTextOutput("verbatRoot", placeholder = FALSE)),
          column(4, verbatimTextOutput("verbatPath1", placeholder = FALSE)),
          column(4, verbatimTextOutput("verbatPath2", placeholder = FALSE))
        ),
        h4("Adjust Visualization"),
        fluidRow(
          column(3, sliderInput("trPlotNodeSize", "Node Size", min = 0.1, max = 3, value = 1)),
          column(3, sliderInput("trPlotCellSize", "Cell Size", min = 0.1, max = 2, value = 1)),
          column(3, sliderInput("trPlotCellAlpha", "Cell Alpha", min = 0.1, max = 1, value = 0.5)),
          column(3, sliderInput("trPlotSegSize", "Trajectory Width", min = 0.1, max = 1, value = 0.5))
        ),
        fluidRow(
          column(3, sliderInput("trPlotNodeText", "Node Text Size", min = 0.5, max = 3, value = 1)),
          column(3, sliderInput("trPlotCellStroke", "Cell Stroke", min = 0.1, max = 2, value = 0.7)),
          column(3, radioButtons("trPlotCellColor", "Color Cells",
            choices = list(
              "Pseudotime" = "pTime",
              "Cell Annotations" = "anno"
            ), selected = "pTime"
          ))
        )
      )
    )
  )


  if (inputType == "Monocle3") {
    # Server
    server <- function(input, output, session) {
      selected_paths <- reactiveVal(list())
      path1_label <- reactive({ input$path1_label })
      path2_label <- reactive({ input$path2_label })
      root_node <- reactiveVal(NULL)
      path1 <- reactiveVal(NULL)
      path2 <- reactiveVal(NULL)
      showSubTrPlot <- reactiveVal(FALSE)
      nodepTime <- reactiveVal(NULL)

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

        if (input$trPlotCellColor == "pTime") {
          colName <- pseudotime_colname
          legendTitle <- "Inferred Pseudotime"
          annotation_data[[colName]] <- round(annotation_data[[colName]], 3)
        } else if (input$trPlotCellColor == "anno") {
          colName <- "anno"
          legendTitle <- "Supplied Annotation"
        }

        # Check for infinite pseudotime
        helptext <- ""
        if (any(is.infinite(annotation_data[[pseudotime_colname]]))) {
          helptext <- "Inf Pseudotime detected (Grey)"
          showNotification("This dataset contains 'Inf' pseudotime. Cells are coloured in Grey", duration = 10, type = "warning")
        }

        if (all(c("path1_high", "path2_high", "root_high") %in% colnames(label_coords))) {
          trajectory.map <- ggplot() +
            geom_point(
              data = annotation_data, aes(x = .data$x, y = .data$y, color = .data[[colName]]), size = input$trPlotCellSize, alpha = input$trPlotCellAlpha,
              stroke = input$trPlotCellStroke
            ) +
            geom_segment(data = trajectory_data, aes(x = .data$x_from, y = .data$y_from, xend = .data$x_to, yend = .data$y_to), linewidth = input$trPlotSegSize) +
            geom_point(
              data = label_coords, aes(x = .data$x, y = .data$y), size = input$trPlotNodeSize,
              color = ifelse(label_coords$root_high == "Yes", "red",
                ifelse(label_coords$path1_high == "Yes" & label_coords$path2_high == "Yes", "purple",
                  ifelse(label_coords$path1_high == "Yes", "green",
                    ifelse(label_coords$path2_high == "Yes", "blue", "black")
                  )
                )
              )
            ) +
            geom_text(data = label_coords, aes(x = .data$x, y = .data$y, label = .data$node), vjust = 1.5, hjust = 0.5, size = input$trPlotNodeText) +
            theme_minimal() +
            theme(legend.position = "none") +
            labs(x = "UMAP-1", y = "UMAP-2")

          if (input$trPlotCellColor == "pTime") {
              trajectory.map <- trajectory.map + scale_color_gradientn(colors = viridis_colors)
        } else if (input$trPlotCellColor == "anno") {
            trajectory.map <- trajectory.map + scale_color_manual(values = cell_colors)
        }
              
          ggplotly(trajectory.map, source = "trajectoryPlot", tooltip = colName) %>%
            plotly::layout(dragmode = "lasso") %>%
            config(displaylogo = FALSE)
        }
      })


      # Set Root Node
      observeEvent(input$rootNode, {
        # Get Selected Data
        lasso_data <- event_data("plotly_selected", source = "trajectoryPlot")
        lasso_data <- lasso_data[lasso_data$curveNumber == max(lasso_data$curveNumber), , drop = FALSE]

        # Get the selected point (Unique from layer-1)
        pointNumber <- as.numeric(unique(lasso_data$pointNumber)) + 1

        # Attach the prefix for Monocle3
        root_vector <- paste("Y", pointNumber, sep = "_")

        # Validate
        if (!is.null(lasso_data) && length(root_vector) == 1) {
          root_node(root_vector)

          # Set log message
          showNotification(paste("Selected Root Node:", root_node()), type = "message")
        } else {
          showNotification(paste("Please select only one node as root node."), type = "error")
        }
      })

      output$verbatRoot <- renderText({
        root_node_value <- root_node()

        return(paste0("Selected Root Nodes: ", root_node_value))
      })
      output$verbatPath1 <- renderText({
        path1_highlight <- path1()

        return(paste0("Selected Path1 Nodes:", paste(path1_highlight, collapse = ",")))
      })
      output$verbatPath2 <- renderText({
        path2_highlight <- path2()
        return(paste0("Selected Path2 Nodes:", paste(path2_highlight, collapse = ",")))
      })

      # Set Path-1
      observeEvent(input$Path1, {
        # Get the lasso data
        lasso_data <- event_data("plotly_selected", source = "trajectoryPlot")
        lasso_data <- lasso_data[lasso_data$curveNumber == max(lasso_data$curveNumber), , drop = FALSE]

        # Get the y nodes
        pointNumber <- as.numeric(unique(lasso_data$pointNumber)) + 1

        # Attach Prefix
        path1_vector <- paste("Y", pointNumber, sep = "_")

        # Validate
        if (!is.null(lasso_data) && length(path1_vector) > 2) {
          if (is.null(root_node())) {
            showNotification(paste("Please Select root node first."), type = "warning")
          } else {
            # Check if exist in root
            if (root_node() %in% path1_vector) {
              # return
              path1(path1_vector)

              # Log
              showNotification(paste("Selected nodes for Path1:", paste(path1(), collapse = ", ")), type = "message")
            } else {
              showNotification(paste("Root node is not a part of selected Path1", root_node(), "should exist in Path1"), type = "error")
            }
          }
        } else {
          showNotification(paste("Please Select more than two nodes to make a path."), type = "warning")
        }
      })
      # Set Path-2
      observeEvent(input$Path2, {
        # Get the lasso data
        lasso_data <- event_data("plotly_selected", source = "trajectoryPlot")
        lasso_data <- lasso_data[lasso_data$curveNumber == max(lasso_data$curveNumber), , drop = FALSE]

        # Get the y nodes
        pointNumber <- as.numeric(unique(lasso_data$pointNumber)) + 1

        # Attach Prefix
        path2_vector <- paste("Y", pointNumber, sep = "_")

        # Validate
        if (!is.null(lasso_data) && length(path2_vector) > 2) {
          if (is.null(root_node())) {
            showNotification(paste("Please Select root node first."), type = "warning")
          } else {
            # Check if exist in root
            if (root_node() %in% path2_vector) {
              if (length(path2_vector[path1() %in% path2_vector]) <= 1) {
                path2(path2_vector)
                # Log
                showNotification(paste("Selected nodes for Path2:", paste(path2(), collapse = ", ")), type = "message")
              } else {
                showNotification(paste("Path2 cannot share nodes with Path1, other than root node", root_node()), type = "error")
              }
            } else {
              showNotification(paste("Root node is not a part of selected Path2", root_node(), "should exist in Path2"), type = "error")
            }
          }
        } else {
          showNotification(paste("Please Select more than two nodes to make a path."), type = "warning")
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

        if (input$trPlotCellColor == "pTime") {
          colName <- pseudotime_colname
        } else if (input$trPlotCellColor == "anno") {
          colName <- "anno"
        }

        sub.trajectory.map <- ggplot() +
          geom_point(
            data = annotation_data_sub, aes(x = .data$x, y = .data$y, color = .data[[colName]]), size = input$trPlotCellSize, alpha = input$trPlotCellAlpha,
            stroke = input$trPlotCellStroke
          ) +
          geom_segment(data = trajectory_data_sub, aes(x = .data$x_from, y = .data$y_from, xend = .data$x_to, yend = .data$y_to), linewidth = input$trPlotSegSize) +
          geom_point(
            data = label_coords_sub, aes(x = .data$x, y = .data$y), size = input$trPlotNodeSize,
            color = ifelse(label_coords_sub$root_high == "Yes", "red",
              ifelse(label_coords_sub$path1_high == "Yes" & label_coords_sub$path2_high == "Yes", "purple",
                ifelse(label_coords_sub$path1_high == "Yes", "green",
                  ifelse(label_coords_sub$path2_high == "Yes", "blue", "black")
                )
              )
            )
          ) +
          geom_text(data = label_coords_sub, aes(x = .data$x, y = .data$y, label = .data$node), vjust = 1.5, hjust = 0.5, size = input$trPlotNodeText) +
          theme_minimal() +
          xlab("UMAP-1") +
          ylab("UMAP-2")
        
        if (input$trPlotCellColor == "pTime") {
            sub.trajectory.map <- sub.trajectory.map + scale_color_gradientn(colors = viridis_colors)
      } else if (input$trPlotCellColor == "anno") {
          sub.trajectory.map <- sub.trajectory.map + scale_color_manual(values = cell_colors)
      }

        ggplotly(sub.trajectory.map) %>% config(displaylogo = FALSE)
      })
      
      # Save button behavior
      observeEvent(input$save, {
        # Get variables
        if (length(path1()) > 2 & length(path2()) > 2 & length(root_node()) == 1) {
          # Set Reactive
          selected_paths(list(root = root_node(), path1 = path1(), path2 = path2(),
                              path1_label = path1_label(), path2_label = path2_label()))

          # Check before close
          if (length(selected_paths()) != 5) {
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
          if (length(selected_paths()) != 5) {
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
  } else {
    stop("Currently only monocle3 is supported")
  }

  selection <- suppressMessages(runApp(shinyApp(ui, server)))

  if (is.null(selection)) {
    return(NULL)
  } else {
    return(selection)
  }
}
