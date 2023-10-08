# library(shiny)
# library(ggplot2)
# library(UpSetR)
# library(DT)
# 
# # Define the UI
# ui <- fluidPage(
#     titlePanel("Shiny Visualization for scmp.obj"),
#     
#     fluidRow(
#         column(4,
#                selectInput("paths", "Choose a variable", colnames(scmp.obj@siggenes@summary)),
#                selectInput("gene", "Choose a variable 2", choices = NULL)
#         ),
#         column(4,
#                plotlyOutput("dynamicPlot")
#         ),
#         column(4,
#                plotOutput("boxi")
#                )
#     ),
#     
#     fluidRow(
#         column(6,
#                dataTableOutput("mytable")
#         ),
#         column(6,
#                plotOutput("dynamicPlot2"))
#     )
# )
# 
# 
# # Define server logic
# server <- function(input, output, session) {
#     
#     gene.counts <- reactiveVal()
#     observe({
#         updateSelectizeInput(session, "gene", choices = scmp.obj@siggenes@summary[, input$paths], server = TRUE)
#     })
#     
#     observe({
#         blk.counts <- scmp.obj@scTFit@dat
#         filtered_data <- blk.counts[rownames(blk.counts) == input$gene,]
#         gene.counts(filtered_data)  # Correctly updating reactiveVal
#     })
#     
#     output$dynamicPlot <- renderPlotly({
#         req(gene.counts())  # Ensuring it has a valid value
#         ggplot_obj <- ggplotly(sc.PlotGroups(data = gene.counts(), dis = scmp.obj@scTFit@dis,
#                                              edesign =  scmp.obj@scTFit@edesign,
#                                              groups.vector = scmp.obj@scTFit@groups.vector, main = input$gene))
#         config(ggplot_obj, displaylogo = FALSE)
#     })
#     
#     output$dynamicPlot2 <- renderPlot({
#         
#         gene.list <- lapply(colnames(scmp.obj@siggenes@summary), function(path){
#             return(na.omit(scmp.obj@siggenes@summary[[path]]))
#         })
#         names(gene.list) <- colnames(scmp.obj@siggenes@summary)
#         
#         upset(fromList(gene.list), main.bar.color = "#F58A53",
#               matrix.color = "#15918A",
#               line.size = 1.5, point.size = 3,
#               sets.x.label = "Number of Features", sets.bar.color = "#EE446F")
#     })
#     
#     output$mytable <- renderDT({
#         datatable(
#             as.data.frame(round(showSol(scmp.obj, view = F, return = T), 3)),
#             options = list(
#                 scrollX = TRUE,
#                 scrollY = "400px",
#                 scrollCollapse = TRUE,
#                 rowCallback = JS(
#                     "function(row, data) {",
#                     sprintf("  if(data[0] == '%s'){", input$gene), 
#                     "    $(row).css('background-color', '#FFDD00');",
#                     "  }",
#                     "}"
#                 )
#             )
#         )
#     
#     })
#     
#     output$boxi <- renderPlot({
#         
#         req(gene.counts()) 
#         
#         
#         boxplot(gene.counts())
#         
# 
#     })
#     
# }
# 
# # Run the app
# shinyApp(ui, server)
# 
