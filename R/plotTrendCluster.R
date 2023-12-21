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
#' @param hclust_agglo_method description
#' @param distance description
#' @param smoothness description
#' @param logs Whether to plot log of counts
#' @param logType Log type required
#' @param includeInflu description
#' @param k description
#' @param result description
#' @param significant description
#' @param parallel description
#' @param verbose description
#' @param summary_mode description
#' @param kmeans_iter_max description
#' @param mclust_k description
#' @param pseudoCount description
#' @param scale_counts description
#' @param summary_mode description
#'
#' @import ggplot2
#' @importFrom stats complete.cases cutree hclust
#' @importFrom RColorConesa getConesaColors
#' @importFrom mclust Mclust
#' @importFrom stats as.dist cor kmeans
#' @return Generates a plot.
#' @export
plotTrendCluster <- function(scmpObj, geneSet, xlab = "Pooled Pseudotime", ylab = "Pseudobulk Expression",
                             cluster_by = "coeff",
                             cluster_method = "hclust",
                             hclust_agglo_method = "ward.D",
                             kmeans_iter_max = 500,
                             distance = "cor",
                             summary_mode = "median",
                             logs = TRUE, logType = "log",
                             smoothness = 0.01, k = 9, includeInflu = TRUE,
                             result = "plot",
                             mclust_k = FALSE,
                             verbose = FALSE,
                             pseudoCount = 1,
                             scale_counts = TRUE,
                             significant = FALSE, parallel = FALSE) {
  # Global vars
  scmp_clusters <- "scmp_clusters"
  pooled.time <- "pooled.time"
  pb.counts <- "pb.counts"
  path <- "path"

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
  assert_that(any(cluster_method %in% c("hclust", "kmeans", "Mclust")),
    msg = paste(
      paste0("'", cluster_method, "'"), "is not a valid method. Please use one of",
      paste(c("hclust", "kmeans", "Mclust"), collapse = ", ")
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

  # Remove NA
  time <- scmpObj@design@alloc[, scmpObj@param@bin_pseudotime_colname, drop = TRUE]
  repvect <- scmpObj@design@alloc[, "Replicate", drop = TRUE]
  dat <- as.data.frame(cluster_matrix_input[, (ncol(cluster_matrix_input) - length(time) + 1):ncol(cluster_matrix_input)])
  count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
  clusterdata <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= length(unique(repvect))), ]

  # NA Treatment
  if (any(is.na(clusterdata))) {
    if (cluster_method == "kmeans" || cluster_method == "Mclust") {
      clusterdata[is.na(clusterdata)] <- 0
      mean.replic <- function(x) {
        tapply(as.numeric(x), repvect, mean, na.rm = TRUE)
        MR <- t(apply(clusterdata, 1, mean.replic))
        # por si acaso todos los valores de una misma r?plica son NAs:
        if (any(is.na(MR))) {
          row.mean <- t(apply(MR, 1, mean, na.rm = TRUE))
          MRR <- matrix(row.mean, nrow(MR), ncol(MR))
          MR[is.na(MR)] <- MRR[is.na(MR)]
        }
        data.noNA <- matrix(NA, nrow(clusterdata), ncol(clusterdata))
        u.repvect <- unique(repvect)
        for (i in 1:nrow(clusterdata)) {
          for (j in 1:length(u.repvect)) {
            data.noNA[i, repvect == u.repvect[j]] <- MR[i, u.repvect[j]]
          }
        }
        clusterdata <- data.noNA
      }
    }
  }

  # Clustering
  if (!is.null(clusterdata)) {
    k <- min(k, nrow(dat), na.rm = TRUE)

    if (cluster_method == "hclust") {
      if (distance == "cor") {
        dcorrel <- matrix(rep(1, nrow(clusterdata)^2), nrow(clusterdata), nrow(clusterdata)) - cor(t(clusterdata),
          use = "pairwise.complete.obs"
        )
        clust <- hclust(as.dist(dcorrel), method = hclust_agglo_method)
        c.algo.used <- paste(cluster_method, "cor", hclust_agglo_method, sep = "_")
      } else {
        clust <- hclust(dist(clusterdata, method = distance), method = hclust_agglo_method)
        c.algo.used <- paste(cluster_method, distance,
          hclust_agglo_method,
          sep = "_"
        )
      }
      cut <- cutree(clust, k = k)
    } else if (cluster_method == "kmeans") {
      cut <- kmeans(clusterdata, k, kmeans_iter_max)$cluster
      c.algo.used <- paste("kmeans", k, kmeans_iter_max, sep = "_")
    } else if (cluster_method == "Mclust") {
      if (mclust_k) {
        my.mclust <- Mclust(clusterdata, verbose = FALSE)
        k <- my.mclust$G
      } else {
        my.mclust <- Mclust(clusterdata, k, verbose = FALSE)
      }
      cut <- my.mclust$class
      c.algo.used <- paste("Mclust", k, sep = "_")
    }
  } else {
    return(
      warning("Impossible to compute hierarchical clustering")
    )
  }
  
  if()

  # Convert to df
  clusters_df <- as.data.frame(cut)
  colnames(clusters_df) <- scmp_clusters
  clusters_df[["feature_id"]] <- rownames(clusters_df)

  # If paralle requested
  if (parallel) {
    os_name <- get_os()
    if (os_name == "windows") {
      numCores <- 1
      warning("Currently, we only support sequential processing on windows based systems...")
    } else {
      numCores <- availableCores() - 1
    }
    if (verbose) {
      message(paste("Running with", numCores, "cores..."))
    }
  } else {
    numCores <- 1
  }

  curve.line.point.list <- mclapply(unique(clusters_df[["feature_id"]]), function(gene_i) {
    # Plot data
    plt <- plotTrend(
      scmpObj = scmpObj,
      feature_id = gene_i,
      smoothness = smoothness,
      xlab = xlab, ylab = ylab,
      logs = logs, logType = logType,
      pseudoCount = pseudoCount,
      significant = significant,
      summary_mode = summary_mode
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
  }, mc.cores = numCores)

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
  point.df[[scmp_clusters]] <- as.factor(point.df[[scmp_clusters]])
  line.df[[scmp_clusters]] <- as.factor(line.df[[scmp_clusters]])
  curve.df[[scmp_clusters]] <- as.factor(curve.df[[scmp_clusters]])
  point.df[["feature_id"]] <- as.factor(point.df[["feature_id"]])
  line.df[["feature_id"]] <- as.factor(line.df[["feature_id"]])
  curve.df[["feature_id"]] <- as.factor(curve.df[["feature_id"]])

  # Set names
  colnames(point.df) <- c("feature_id", pooled.time, pb.counts, path, scmp_clusters)
  colnames(line.df) <- c("feature_id", pooled.time, path, pb.counts, scmp_clusters)
  colnames(curve.df) <- c("feature_id", pooled.time, pb.counts, path, scmp_clusters)

  # Drop features
  line.df <- line.df[, colnames(line.df) != "feature_id", drop = FALSE]
  curve.df <- curve.df[, colnames(curve.df) != "feature_id", drop = FALSE]

  # Take medians
  line.df <- line.df %>%
    group_by(!!sym(path), !!sym(scmp_clusters), !!sym(pooled.time)) %>%
    summarize(pb.counts.mean = median(!!sym(pb.counts)), .groups = "drop")
  curve.df <- curve.df %>%
    group_by(!!sym(path), !!sym(scmp_clusters), !!(pooled.time)) %>%
    summarize(pb.counts.mean = median(!!sym(pb.counts)), .groups = "drop")

  if (result == "plot") {
    # Plot
    p <- ggplot() +
      geom_point(
        data = point.df, aes(x = .data$pooled.time, y = .data$pb.counts, color = .data$path),
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
          x = .data$pooled.time, y = .data$pb.counts,
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
      labs(title = "Gene Expression over Pseudotime", color = "Branching Path") +
      xlab(xlab) +
      ylab(ylab)
    return(p)
  }
}
