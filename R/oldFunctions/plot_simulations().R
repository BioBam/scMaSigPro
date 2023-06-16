# Function to plot simulated data from splatter in lower dimension
plot_simulations <- function(sim.sce,
                             title.2d = NULL,
                             title.3d = NULL,
                             plot2d = FALSE, plot3d = TRUE,
                             colorGroupDiscrete = "Group",
                             colorGroupContinuous = "Step",
                             merge.by = "Cell", assay_type = "rawCounts") {
  # Required Libraries
  suppressPackageStartupMessages(require(scran))
  suppressPackageStartupMessages(require(scater))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(plotly))

  # Calculate Log-normalized Counts
  sim.sce <- logNormCounts(sim.sce,
    transform = "log",
    assay.type = assay_type
  )

  # Compute PCA
  sim.sce <- runPCA(sim.sce)

  # Extract Components for 3D visulization
  threeD_PCA <- as.data.frame(reducedDim(sim.sce)[, c(1:3)])

  # Add cells to table
  threeD_PCA[[merge.by]] <- rownames(threeD_PCA)

  # Merge with Metadata
  threeD_PCA <- merge(threeD_PCA, as.data.frame(colData(sim.sce)), by = merge.by)

  # Require Plotly
  if (plot3d == T) {
    p1 <- plot_ly(threeD_PCA,
      x = ~PC1, y = ~PC2, z = ~PC3, color = ~ .data[[colorGroupContinuous]], alpha = 0.7, mode = "markers",
      size = 2
    ) %>%
      layout(title = paste(title.3d, "color by", colorGroupContinuous), plot_bgcolor = "#e5ecf6")

    p2 <- plot_ly(threeD_PCA,
      x = ~PC1, y = ~PC2, z = ~PC3, color = ~ .data[[colorGroupDiscrete]], mode = "markers", alpha = 0.7,
      size = 2
    ) %>% layout(
      title = paste(title.3d, "color by", colorGroupDiscrete), plot_bgcolor = "#e5ecf6",
      legend = list(title = list(text = colorGroupDiscrete))
    )

    print(p1)
    print(p2)
  }

  # Basic PCA Plot
  if (plot2d == T) {
    p3 <- plotPCA(sim.sce, colour_by = colorGroupDiscrete) + ggtitle(paste(title.2d, "color by", colorGroupDiscrete))
    p4 <- plotPCA(sim.sce, colour_by = colorGroupContinuous) + ggtitle(paste(title.2d, "color by", colorGroupContinuous))
    print(ggpubr::ggarrange(p3, p4))
  }
}


plot_cell_time_distribution <- function(sce_object,
                                        cell_col = "Cell", time_col = "Step",
                                        path_col = "Group", path.vec = c("Path1", "Path2")) {
  cell.meta <- as.data.frame(colData(sce_object))

  # Select Columns
  plt.table <- cell.meta[, c(cell_col, time_col, path_col)]

  # Group by
  # plt.table <- plt.table %>% group_by(!!!time_col, !!!path_col) %>% summarise(cluster.members=paste0(!!!cell_col, collapse='|'))
  plt.table <- plt.table %>%
    group_by(.dots = c(time_col, path_col)) %>%
    summarise(cluster.members = paste0(Cell, collapse = "|"))

  # Count number of cells
  plt.table$Num <- apply(plt.table, 1, calc_bin_size)

  # Select Columns
  plt.table <- plt.table[, !colnames(plt.table) %in% "cluster.members"]

  # Separation of plots
  plt.table.1 <- plt.table[plt.table[[path_col]] == path.vec[1], ]
  plt.table.2 <- plt.table[plt.table[[path_col]] == path.vec[2], ]

  # Plot Data
  p <- ggplot(plt.table.1, aes(x = Num)) +
    geom_histogram(
      binwidth = 0.5, ,
      color = "#f68a53", fill = "#f68a53", alpha = 05
    ) +
    geom_vline(aes(xintercept = mean(Num)), linetype = "dashed", color = "#139289") +
    theme_classic() +
    theme(
      legend.position = "none", strip.text = element_text(size = rel(2)),
      axis.text = element_text(size = rel(1)),
      panel.grid.major = element_line(size = 0.7, linetype = "dotted"),
      panel.grid.minor = element_line(size = 0.2)
    ) +
    ggtitle("Distribution of cells per Time-point", path.vec[1]) +
    ylab("Number of cells") +
    xlab("Number of cells associations")

  p2 <- ggplot(plt.table.2, aes(x = Num)) +
    geom_histogram(
      binwidth = 0.5, ,
      color = "#f68a53", fill = "#f68a53", alpha = 05
    ) +
    geom_vline(aes(xintercept = mean(Num)), linetype = "dashed", color = "#139289") +
    theme_classic() +
    theme(
      legend.position = "none", strip.text = element_text(size = rel(2)),
      axis.text = element_text(size = rel(1)),
      panel.grid.major = element_line(size = 0.7, linetype = "dotted"),
      panel.grid.minor = element_line(size = 0.2)
    ) +
    ggtitle("Distribution of cells per Time-point", path.vec[2]) +
    ylab("Number of cells") +
    xlab("Number of cells associations")

  print(p)
  print(p2)
}
