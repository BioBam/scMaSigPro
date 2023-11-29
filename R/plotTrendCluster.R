#' Plot Groups Function
#'
#' This function generates plots based on various parameters.
#'
#' @param scmpObj object of class scmpObj
#' @param groupBy A vector of maximum length 6 or a single feature of ID
#' @param xlab X-axis label. Default is "Pooled Pseudotime".
#' @param ylab Y-axis label. Default is "Pseudobulk Expression".
#' @param smoothness description
#' @param logs Whether to plot log of counts
#' @param logType Log type required
#' @param includeInflu description
#'
#' @import ggplot2
#' @importFrom RColorConesa getConesaColors
#' @return Generates a plot.
#' @export

plotTrendCluster <-
  function(scmpObj,
           xlab = "Pooled Pseudotime",
           ylab = "Pseudobulk Expression",
           groupBy = "coeff",
           logs = TRUE,
           logType = "log",
           smoothness = 0.01,
           includeInflu = TRUE) {
    # Get the names of the genes that are significant
    sig_gene_list <- scmpObj@sig.genes@sig.genes
    sig_genes <- unlist(sig_gene_list)
    cluster_list <- scmpObj@sig.genes@feature.clusters

    # Create Data list
    trend.data.list <- lapply(sig_genes, function(gene_i, group_by = groupBy, sigGeneList = sig_gene_list,
                                                  clusterList = cluster_list) {
      # Create the plot
      plt <- plotTrend(
        scmpObj = scmpObj,
        feature_id = gene_i,
        smoothness = smoothness,
        xlab = xlab, ylab = ylab,
        logs = logs, logType = logType
      )

      # Extract point data
      point.data <- plt$layers[[1]][["data"]]

      # Extract
      if (group_by == "feature") {
        trend.data <- plt$layers[[2]][["data"]]
        trend.data[["feature"]] <- gene_i
        clusterList <- clusterList[["sigCounts"]]
        clusterLabel <- clusterList[[gene_i]]
        trend.data[["scmpCluster"]] <- paste("cluster", clusterLabel)
      } else if (group_by == "coeff") {
        trend.data <- plt$layers[[3]][["data"]]
        trend.data[["feature"]] <- gene_i
        clusterList <- clusterList[["sigCoeff"]]
        clusterLabel <- clusterList[[gene_i]]
        trend.data[["scmpCluster"]] <- paste("cluster", clusterLabel)
      }

      # Reset columns
      colnames(trend.data) <- c("x_axis", "y_axis", "path", "feature_id", "scmpCluster_id")

      return(trend.data)
    })

    # Get list
    trend.data <- do.call("rbind", trend.data.list)
    rownames(trend.data) <- NULL

    # Assuming feature_id and cluster_id are factors
    p <- ggplot(
      data = trend.data,
      aes(
        x = .data$x_axis, y = .data$y_axis,
        group = interaction(.data$feature_id, .data$path),
        color = .data$path, shape = .data$path, linetype = .data$path
      )
    ) +
      geom_line(
        linewidth = 0.4
      ) + # Draw lines
      geom_point(
        size = 1, alpha = 0.5, stroke = 1
      ) + # Draw points
      facet_wrap(~ .data$scmpCluster_id, scales = "free_y") + # Create a panel for each cluster_id
      scale_color_manual(values = colorConesa(length(unique(trend.data$path)))) + # Custom colors for paths
      theme_classic(base_size = 10) +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, angle = 0),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis text if necessary
      ) +
      labs(title = "Gene Expression over Pseudotime", color = "Path") +
      xlab(xlab) +
      ylab(ylab)
    # guides(color = guide_legend(title = "Path")) # Adjust legend for paths
    return(p)
  }
