#' Plot Groups Function
#'
#' This function generates plots based on various parameters.
#'
#' @param scmpObj object of class scmpObj
#' @param xlab X-axis label. Default is "Pooled Pseudotime".
#' @param ylab Y-axis label. Default is "Pseudobulk Expression".
#' @param cluster_by description
#' @param geneSet description
#' @param cluster_method description
#' @param hclust.agglo_method description
#' @param distance description
#' @param smoothness description
#' @param logs Whether to plot log of counts
#' @param logType Log type required
#' @param includeInflu description
#' @param k description
#' @param result description
#' @param significant description
#'
#' @import ggplot2
#' @importFrom stats complete.cases cutree hclust
#' @importFrom RColorConesa getConesaColors
#' @return Generates a plot.
#' @export
plotTrendCluster <- function(scmpObj, geneSet, xlab = "Pooled Pseudotime", ylab = "Pseudobulk Expression",
                             cluster_by = "coeff",
                             cluster_method = "hclust", logs = TRUE, logType = "log",
                             smoothness = 0.01, k = 9, includeInflu = TRUE, distance = "cor",
                             hclust.agglo_method = "ward.D",
                             result = "plot",
                             significant = FALSE) {
  # Check if the gene set exists
  assert_that(any(geneSet %in% c(names(scmpObj@sig.genes@sig.genes), "intersect", "union")),
    msg = paste(
      paste0("'", geneSet, "'"), "does not exist. Please use one of",
      paste(c(names(scmpObj@sig.genes@sig.genes), "intersect", "union"), collapse = ", ")
    )
  )
  assert_that(any(cluster_by %in% c("coeff", "counts")),
    msg = paste(
      paste0("'", cluster_by, "'"), "is not a valid option. Please use one of",
      paste(c("coeff", "counts"), collapse = ", ")
    )
  )
  assert_that(any(cluster_method %in% c("hclust", "kmeans")),
    msg = paste(
      paste0("'", cluster_method, "'"), "is not a valid method. Please use one of",
      paste(c("hclust", "kmeans"), collapse = ", ")
    )
  )
  assert_that(any(result %in% c("return", "plot")),
    msg = paste(
      paste0("'", result, "'"), "is not a valid option. Please use one of",
      paste(c("return", "plot"), collapse = ", ")
    )
  )

  # Get gene set vector
  if (geneSet == "intersect") {
    gene_set_vector <- Reduce(intersect, scmpObj@sig.genes@sig.genes)
  } else if (geneSet == "union") {
    gene_set_vector <- Reduce(union, scmpObj@sig.genes@sig.genes)
  } else {
    gene_set_vector <- scmpObj@sig.genes@sig.genes[[geneSet]]
  }

  # Extract data based on 'cluster_by'
  if (cluster_by == "counts") {
    # Extract bulk counts
    cluster_matrix_input <- as.matrix(showSigProf(scmpObj, includeInflu = includeInflu))
    cluster_matrix_input <- cluster_matrix_input[rownames(cluster_matrix_input) %in% gene_set_vector, , drop = FALSE]
  } else if (cluster_by == "coeff") {
    # Extract coefficients
    cluster_matrix_input <- as.matrix(showCoeff(scmpObj, includeInflu = includeInflu))
    cluster_matrix_input <- cluster_matrix_input[rownames(cluster_matrix_input) %in% gene_set_vector, , drop = FALSE]
  }

  # Check method
  if (cluster_method == "hclust") {
    cluster_matrix_input[is.na(cluster_matrix_input)] <- 0
    # cluster_matrix_input <- cluster_matrix_input[complete.cases(cluster_matrix_input), , drop = FALSE]

    # Compute the distance matrix
    dist_matrix <- dist(cluster_matrix_input)

    # Perform hierarchical clustering
    hc <- hclust(dist_matrix)

    # Cut the tree to form nine clusters
    clusters <- cutree(hc, k = k)

    # Convert to matrix for mathcing
    clusters_df <- data.frame(
      scmp_clusters = clusters,
      feature_id = names(clusters)
    )

    if (result == "return") {
      cluster.list <- list(clusters)
      names(cluster.list) <- geneSet
      scmpObj@sig.genes@feature.clusters <- cluster.list
      return(scmpObj)
    }

    clusters_df[["scmp_clusters"]] <- paste("Cluster:", clusters_df[["scmp_clusters"]])
    rownames(clusters_df) <- NULL
  } else {
    stop("not coded")
    # else if (cluster_method == "kmeans") {
    #     # Standardizing the data can be important for k-means
    #     standardized_data <- scale(sig.element)
    #
    #     # Run k-means clustering
    #     clustering.result <- kmeans(standardized_data, centers = numClus, nstart = 25)
    #     # nstart parameter is used to set the number of random sets chosen
    #
    #     # Extract Vector
    #     cluster.vector <- as.vector(clustering.result$cluster)
    #     names(cluster.vector) <- rownames(sig.element)
    # } else if (cluster_method == "mclust") {
    #     # Run Mclust
    #     # Mclust automatically determines the optimal number of clusters
    #     # You can specify the number of clusters if desired using the G parameter
    #     clustering.result <- Mclust(sig.element, G = numClus)
    #
    #     # Extract the best model and get cluster assignments
    #     bestModel <- summary(clustering.result, parameters = TRUE)
    #     cluster.vector <- as.vector(as.integer(bestModel$classification))
    #
    #     # Optionally, you can extract the BIC values and other model details
    #     # bic_values <- clustering.result$BIC
    #
    #     # Assign row names
    #     names(cluster.vector) <- rownames(sig.element)
    # }
  }

  curve.line.point.list <- mclapply(unique(clusters_df[["feature_id"]]), function(gene_i) {
    # Plot data
    plt <- plotTrend(
      scmpObj = scmpObj,
      feature_id = gene_i,
      smoothness = smoothness,
      xlab = xlab, ylab = ylab,
      logs = logs, logType = logType,
      significant = significant
    )

    # Extract layers
    point.data <- plt$layers[[1]][["data"]]
    line.data <- plt$layers[[2]][["data"]]
    curve.data <- plt$layers[[3]][["data"]]

    # Add gene names
    point.data[["feature_id"]] <- gene_i
    line.data[["feature_id"]] <- gene_i
    curve.data[["feature_id"]] <- gene_i

    # Add cluster levelinfo
    line.data <- merge(line.data, clusters_df, by = "feature_id")
    point.data <- merge(point.data, clusters_df, by = "feature_id")
    curve.data <- merge(curve.data, clusters_df, by = "feature_id")

    return(
      list(
        curve.data = curve.data,
        point.data = point.data,
        line.data = line.data
      )
    )
  }, mc.cores = availableCores())

  # Sepatarte List
  line.list <- lapply(curve.line.point.list, function(x) {
    return(x[[3]])
  })
  point.list <- lapply(curve.line.point.list, function(x) {
    return(x[[2]])
  })
  curve.list <- lapply(curve.line.point.list, function(x) {
    return(x[[1]])
  })

  # Bind
  point.df <- do.call("rbind", point.list)
  line.df <- do.call("rbind", line.list)
  curve.df <- do.call("rbind", curve.list)

  # Create factors
  point.df[["scmp_clusters"]] <- as.factor(point.df[["scmp_clusters"]])
  line.df[["scmp_clusters"]] <- as.factor(line.df[["scmp_clusters"]])
  curve.df[["scmp_clusters"]] <- as.factor(curve.df[["scmp_clusters"]])
  point.df[["feature_id"]] <- as.factor(point.df[["feature_id"]])
  line.df[["feature_id"]] <- as.factor(line.df[["feature_id"]])
  curve.df[["feature_id"]] <- as.factor(curve.df[["feature_id"]])

  # Set names
  colnames(point.df) <- c("feature_id", "pooled.time", "pb.counts", "path", "scmp_clusters")
  colnames(line.df) <- c("feature_id", "pooled.time", "path", "pb.counts", "scmp_clusters")
  colnames(curve.df) <- c("feature_id", "pooled.time", "pb.counts", "path", "scmp_clusters")

  # Drop features
  line.df <- line.df[, colnames(line.df) != "feature_id", drop = FALSE]
  curve.df <- curve.df[, colnames(curve.df) != "feature_id", drop = FALSE]

  # Take medians
  line.df <- line.df %>%
    group_by(path, scmp_clusters, pooled.time) %>%
    summarize(pb.counts.mean = median(pb.counts), .groups = "drop")
  curve.df <- curve.df %>%
    group_by(path, scmp_clusters, pooled.time) %>%
    summarize(pb.counts.mean = median(pb.counts), .groups = "drop")

  if (result == "plot") {
    # Plot
    p <- ggplot() +
      geom_point(
        data = point.df, aes(x = pooled.time, y = pb.counts, color = path),
        fill = "#102C57", alpha = 0.5, size = 0.5, stroke = 0.5, shape = 21
      ) +
      geom_line(
        data = line.df,
        aes(
          x = .data$pooled.time, y = .data$pb.counts.mean,
          color = .data$path, group = .data$path,
        ), linetype = "solid", linewidth = 0.5
      ) +
      geom_line(
        data = curve.df,
        aes(
          x = .data$pooled.time, y = .data$pb.counts.mean,
          color = .data$path, group = .data$path,
        ), linetype = "dotted", linewidth = 0.5
      ) +
      facet_wrap(~ .data$scmp_clusters, scales = "free_y") + # Create a panel for each cluster_id
      scale_color_manual(values = colorConesa(length(unique(point.df$path)))) + # Custom colors for paths
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
    return(p)
  }
}
