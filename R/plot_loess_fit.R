#' Plot Spline Smoothers
#'
#' @param sce_obj Object of the \link[SingleCellExperiment]{SingleCellExperiment} class.
#' @param gene_name Name of the gene to be plotted.
#' @param log Whether to plot log of the counts.
#' @param assay_name Name of the assay in the \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param time_col Column name having pseudotime information.
#' @param path_col Column name having path information.
#' @param min_exp Minimum expression to plot.
#' @param plt_subtitle Sub-title to add to the plot.
#' @param dfreedom Degrees of Freedom for loess fit.
#' @param span Smoothing parameter for the loess fit.
#

# Plot B-Spline
plot_loess_fit <- function(sce_obj, gene_name, log = F,
                           assay_name = "TrueCounts",
                           time_col = "Step",
                           min_exp = 0,
                           path_col = "Group",
                           plt_subtitle = NULL,
                           dfreedom = 3,
                           span = 0.75) {
  # Test
  if (length(unique(rownames(sce_obj) %in% gene_name)) == 1) {
    stop("gene names does not exist")
  }

  # Extract the colData
  cell.data <- as.data.frame(colData(sce_obj))[, c(time_col, path_col), drop = F]
  cell.data$cell_id <- rownames(cell.data)
  
  # Extract the matrix depending upon the class of the matrix
  t.counts.mtx <- sce_obj@assays@data[[assay_name]]
  
  # Long format
  suppressWarnings(t.counts.tab <- melt(as.matrix(t.counts.mtx)))
  colnames(t.counts.tab) <- c("gene_name", "cell_id", "counts_value")

  # Select the gene of interest
  plt.tab <- t.counts.tab[t.counts.tab$gene_name == gene_name, , drop = F]

  # Merge with time
  plt.tab <- merge(plt.tab, cell.data, by = "cell_id")

  if (log == T) {
    plt.tab$counts_value <- log(plt.tab$counts_value)
  }

  # return(plt.tab)
  plt.tab <- plt.tab[plt.tab$counts_value >= min_exp, ]
  
  # Base Plot
  loess_plot <- ggplot() +
    theme_classic() +
    ggtitle(
      label = paste0("Gene Trend Plot: ", gene_name),
      subtitle = paste0("df: ", dfreedom, " | ", "Span : ", span, " | ", plt_subtitle)
    ) +
    xlab(paste0("Continuum: ", time_col)) +
    ylab(assay_name) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.key.size = unit(1, "cm"),
      axis.text = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.2)),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.3)
    )
  if (log == T) {
    loess_plot <- loess_plot + ylab(paste0("log(", assay_name, ")"))
  }

  # Add layers
  for (i in levels(as.factor(unique(cell.data[[path_col]])))) {
    # Subset
    plt.tab.path <- plt.tab[plt.tab[[path_col]] == i, , drop = F]

    # Add Points
    loess_plot <- loess_plot +
      geom_point(
        data = plt.tab.path,
        aes(x = .data[[time_col]], y = .data[["counts_value"]], color = .data[[path_col]]),
        alpha = 0.75, fill = "#fdc659", size = rel(0.6)
      )

    # Add Smoother
    loess_plot <- loess_plot +
      geom_smooth(
        formula = y ~ x, data = plt.tab.path, linewidth = 0.75,
        alpha = 0.75, aes(
          x = .data[[time_col]],
          y = .data[["counts_value"]],
          color = .data[[path_col]]
        ),
        fill = "#FFFFCC", method = "loess", se = T,
        method.args = list(
          degree = dfreedom, span = span,
          normalize = F,
          family = "symmetric"
        )
      )
  }

  # Additional Formatting
  loess_plot <- loess_plot + scale_color_brewer(palette = "Dark2")
  loess_plot <- loess_plot + scale_x_continuous(
    breaks = round(seq(
      from = min(plt.tab[[time_col]]),
      to = max(plt.tab[[time_col]]), length.out = 10
    ))
  )
  loess_plot <- loess_plot + scale_y_continuous(
    breaks = round(seq(
      from = min(plt.tab$counts_value),
      to = max(plt.tab$counts_value), length.out = 5
    ))
  )
  return(loess_plot)
}
