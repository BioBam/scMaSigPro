#' @title Plot multiple trends of the multiple genes.
#'
#' @description
#' Plot trends of multiple genes (clustered) across the binned pseudotime.
#'
#' @import ggplot2
#' @importFrom stats complete.cases cutree hclust
#' @importFrom mclust Mclust
#' @importFrom stats as.dist cor kmeans
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param xlab X-axis label. (Default is "Pooled Pseudotime")
#' @param ylab Y-axis label. (Default is "Pseudobulk Expression")
#' @param plot Whether to plot 'coeff' or 'counts'. (Default is 'counts')
#' @param smoothness How smooth the trend should be. Setting to
#' higher values will result in more linear trends. (Default is 0.01)
#' @param logs Whether to log transform counts. (Default is TRUE)
#' @param logType How to log transform the values. Available options 'log',
#' 'log2', 'log10'. (Default is 'log')
#' @param includeInflu Include gene only if it has influential data.
#' (Default is TRUE)
#' @param significant Include gene only if the models are significant based on
#' \code{scMaSigPro::sc.filter()}. (Default is TRUE)
#' @param parallel Use forking process to run parallelly. (Default is FALSE)
#' (Currently, Windows is not supported)
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' @param summary_mode Compress the expression values per replicate (if present)
#'  per binned pseudotime point. Default is 'median'. Other option 'mean'
#' @param pseudoCount Add a pseudo-count before taking the log. (Default is 1)
#' @param curves Whether to plot the fitted curves. (Default is TRUE)
#' @param lines Whether to plot the lines. (Default is FALSE)
#' @param points Whether to plot the points. (Default is TRUE)
#'
#' @return ggplot2 plot object.
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @export
plotTrendCluster <- function(scmpObj,
                             xlab = "Pooled Pseudotime",
                             ylab = "log(Pseudobulk Expression)",
                             plot = "counts",
                             summary_mode = "median",
                             logs = TRUE, logType = "log",
                             smoothness = 1,
                             includeInflu = TRUE,
                             verbose = TRUE,
                             pseudoCount = 1,
                             significant = FALSE,
                             curves = TRUE,
                             lines = FALSE,
                             points = TRUE,
                             parallel = FALSE) {
  # # # Debugg
  # scmpObj <- multi_scmp_ob_A
  # xlab <- "Pooled Pseudotime"
  # ylab <- "Pseudobulk Expression"
  # plot <- "counts"
  # summary_mode <- "median"
  # logs <- TRUE
  # logType <- "log"
  # smoothness <- 1
  # includeInflu <- TRUE
  # verbose <- TRUE
  # pseudoCount <- 1
  # significant <- FALSE
  # c <- curves <- TRUE
  # l <- lines <- TRUE
  # p <- points <- TRUE
  # parallel <- FALSE


  # Global vars
  scmp_clusters <- "scmp_clusters"
  feature_id <- "feature_id"
  # Offset
  offset_vector <- scmpObj@Design@offset

  # Check
  assertthat::assert_that(!isEmpty(scmpObj@Significant@clusters),
    msg = "Please run 'sc.cluster.trend', before plotting cluster trends"
  )

  # Check Assertion
  assertthat::assert_that(curves || lines || points, msg = "At least one of 'curves', 'lines', or 'points' must be TRUE.")

  assertthat::assert_that(any(plot %in% c("coeff", "counts")),
    msg = paste(
      paste0("'", plot, "'"), "is not a valid option. Please use one of",
      paste(c("coeff", "counts"), collapse = ", ")
    )
  )

  # Extract Gene_Set
  gene_set_vector <- names(scmpObj@Significant@clusters)

  # Extarct cluster set
  gene_cluster_vector <- paste("cluster", as.vector(unlist(scmpObj@Significant@clusters)),
    sep = "_"
  )

  # Create df
  gene.cluster.df <- data.frame(
    gene = gene_set_vector,
    clus = gene_cluster_vector
  )
  colnames(gene.cluster.df) <- c(feature_id, scmp_clusters)

  # Calc
  freq.table <- as.data.frame(table(gene.cluster.df[[scmp_clusters]]))
  colnames(freq.table) <- c("cluster", "num")

  # Verbose
  if (verbose) {
    message("Computing trends for each gene, hang tight this can take a few minutes...")
  }

  all.plt.list <- lapply(c(gene.cluster.df[["feature_id"]], "got"), function(gene_i,
                                                                             scmp_obj = scmpObj,
                                                                             smooth = smoothness,
                                                                             log = logs, log_type = logType,
                                                                             pCount = pseudoCount,
                                                                             sig = significant,
                                                                             summary = summary_mode,
                                                                             c = curves,
                                                                             l = lines,
                                                                             p = points) {
    # Run per genes
    plt <- tryCatch(
      {
        # Your original plotTrend call
        plt <- plotTrend(
          scmpObj = scmp_obj,
          feature_id = gene_i,
          smoothness = smooth,
          logs = log, logType = log_type,
          pseudoCount = pCount,
          significant = sig,
          summary_mode = summary,
          curves = c,
          lines = l,
          points = p
        )
        return(plt)
      },
      error = function(e) {
        # In case of an error, return NULL
        return(NULL)
      }
    )

    # Return Plots
    return(plt)
  })

  if (verbose) {
    message("Computed trends for each gene, collapsing trends...")
  }

  # Set names
  names(all.plt.list) <- c(gene.cluster.df[["feature_id"]], "got")

  # Remove Any possible NULL
  plt.list <- Filter(Negate(is.null), all.plt.list)

  # Start Traversing
  data.list <- lapply(plt.list, function(gene_i.plot,
                                         c = curves,
                                         l = lines,
                                         p = points) {
    point.data <- NULL
    line.data <- NULL
    curve.data <- NULL

    # Extract layers
    if (c) {
      curve.data <- gene_i.plot$layers[["curves"]][["data"]]
      colnames(curve.data) <- c("x_axis", "y_axis", "group")
    }
    if (l) {
      line.data <- gene_i.plot$layers[["lines"]][["data"]]
      colnames(line.data) <- c("x_axis", "y_axis", "group")
    }
    if (p) {
      point.data <- gene_i.plot$layers[["points"]][["data"]]
      colnames(point.data) <- c("x_axis", "y_axis", "group")
    }

    return(list(
      points = point.data,
      line = line.data,
      curve = curve.data
    ))
  })

  # Create List as per the clusters
  cluster.data.list <- lapply(unique(gene.cluster.df[[scmp_clusters]]), function(clus,
                                                                                 gene_cluster_df = gene.cluster.df,
                                                                                 cluster = scmp_clusters,
                                                                                 feature = feature_id,
                                                                                 data_list = data.list) {
    # Get gene names with same cluster
    gene_set <- gene_cluster_df[gene_cluster_df[[cluster]] == clus, , drop = FALSE][[feature]]

    # Subset the datalist
    return(data_list[gene_set])
  })

  # Reset names
  names(cluster.data.list) <- paste("cluster", seq_len(length.out = length(cluster.data.list)),
    sep = "_"
  )

  # Df list
  collapsed.df <- lapply(cluster.data.list, function(clus_i.list,
                                                     summary = summary_mode,
                                                     c = curves,
                                                     l = lines,
                                                     p = points) {
    point.df <- NULL
    line.df <- NULL
    curve.df <- NULL
    # Define a summarization function

    if (l) {
      # Extract Sub.list
      line.list <- lapply(clus_i.list, function(line) {
        return(line[["line"]])
      })
      line.df <- do.call("rbind", line.list)
    }

    if (p) {
      point.list <- lapply(clus_i.list, function(point) {
        return(point[["points"]])
      })
      point.df <- do.call("rbind", point.list)
    }

    if (c) {
      curve.list <- lapply(clus_i.list, function(curve) {
        return(curve[["curve"]])
      })
      curve.df <- do.call("rbind", curve.list)
    }



    # Grouping and summarizing with .data pronoun
    if (summary == "mean") {
      if (!is.null(line.df)) {
        line.df <- line.df %>%
          dplyr::group_by(.data$x_axis, .data$group) %>%
          dplyr::summarize(y_axis = mean(.data$y_axis, na.rm = TRUE), .groups = "drop")
      }
      if (!is.null(point.df)) {
        point.df <- point.df %>%
          dplyr::group_by(.data$x_axis, .data$group) %>%
          dplyr::summarize(y_axis = mean(.data$y_axis, na.rm = TRUE), .groups = "drop")
      }
      if (!is.null(curve.df)) {
        curve.df <- curve.df %>%
          dplyr::group_by(.data$x_axis, .data$group) %>%
          dplyr::summarize(y_axis = mean(.data$y_axis, na.rm = TRUE), .groups = "drop")
      }
    } else if (summary == "median") {
      if (!is.null(line.df)) {
        line.df <- line.df %>%
          dplyr::group_by(.data$x_axis, .data$group) %>%
          dplyr::summarize(y_axis = median(.data$y_axis, na.rm = TRUE), .groups = "drop")
      }
      if (!is.null(point.df)) {
        point.df <- point.df %>%
          dplyr::group_by(.data$x_axis, .data$group) %>%
          dplyr::summarize(y_axis = median(.data$y_axis, na.rm = TRUE), .groups = "drop")
      }
      if (!is.null(curve.df)) {
        curve.df <- curve.df %>%
          dplyr::group_by(.data$x_axis, .data$group) %>%
          dplyr::summarize(y_axis = median(.data$y_axis, na.rm = TRUE), .groups = "drop")
      }
    }

    return(list(
      points = point.df,
      line = line.df,
      curve = curve.df
    ))
  })

  # Initialize empty lists to store each type of data frame
  lines_list <- list()
  curves_list <- list()
  points_list <- list()

  # Loop over the list of clusters
  for (i in seq_along(collapsed.df)) {
    cluster_name <- paste0("cluster_", i)

    # Extract 'line', 'curve', and 'points' data frames and add the 'cluster' column
    if (lines) {
      line_df <- collapsed.df[[i]]$line
      line_df$cluster <- cluster_name
      lines_list[[i]] <- line_df
    }

    if (curves) {
      curve_df <- collapsed.df[[i]]$curve
      curve_df$cluster <- cluster_name
      curves_list[[i]] <- curve_df
    }
    if (points) {
      points_df <- collapsed.df[[i]]$points
      points_df$cluster <- cluster_name
      points_list[[i]] <- points_df
    }
  }

  # Combine the data frames of the same type
  if (lines) {
    lines_combined <- do.call(rbind, lines_list) %>% as.data.frame()
    lines_combined <- merge(freq.table, lines_combined, by = "cluster")
    lines_combined[["cluster"]] <- paste0(
      paste("Cluster", stringr::str_split_i(string = lines_combined[["cluster"]], pattern = "_", i = 2), sep = ": "),
      " (", lines_combined[["num"]], " Features)"
    )
  }
  if (curves) {
    curves_combined <- do.call(rbind, curves_list) %>% as.data.frame()
    curves_combined <- merge(freq.table, curves_combined, by = "cluster")
    curves_combined[["cluster"]] <- paste0(
      paste("Cluster", stringr::str_split_i(string = curves_combined[["cluster"]], pattern = "_", i = 2), sep = ": "),
      " (", curves_combined[["num"]], " Features)"
    )
  }
  if (points) {
    points_combined <- do.call(rbind, points_list) %>% as.data.frame()
    points_combined <- merge(freq.table, points_combined, by = "cluster")
    # Fix cluster info
    points_combined[["cluster"]] <- paste0(
      paste("Cluster", stringr::str_split_i(string = points_combined[["cluster"]], pattern = "_", i = 2), sep = ": "),
      " (", points_combined[["num"]], " Features)"
    )
  }

  if (verbose) {
    message("Calculating 'loess', hang tight this can take a while depending on the smoothness...")
  }

  if (sum(offset_vector) != 0) {
    if (lines) {
      lines_combined <- points_combined
    }
  }


  # View(points_combined)
  # View(curves_combined)
  # View(lines_combined)

  # Initiate plotting
  p <- ggplot()

  suppressWarnings(expr = {
    if (points) {
      p <- p + geom_point(
        data = points_combined, aes(x = .data$x_axis, y = .data$y_axis, color = .data$group),
        fill = "#102C57", alpha = 0.5, size = 0.5, stroke = 0.5, shape = 21
      )
    }
    if (lines) {
      p <- p + geom_path(
        data = lines_combined,
        aes(
          x = .data$x_axis, y = .data$y_axis, color = .data$group, group = .data$group,
        ), linetype = "dashed", linewidth = 0.5
      )
    }

    if (curves) {
      p <- p + geom_smooth(
        data = curves_combined,
        se = FALSE,
        formula = y ~ x, span = 0.7,
        method = "loess",
        aes(
          x = .data$x_axis, y = .data$y_axis, color = .data$group, group = .data$group,
        ), linetype = "solid", linewidth = 0.5
      )
    }

    p <- p + facet_wrap(~ .data$cluster, scales = "free_y") + # Create a panel for each cluster_id
      scale_color_manual(values = scmp_colors(length(unique(cDense(scmpObj)[[scmpObj@Parameters@path_col]])))) + # Custom colors for paths
      theme_classic(base_size = 10) +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, angle = 0),
        legend.position = "bottom", legend.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis text if necessary
      ) +
      labs(title = "Gene Expression over Pseudotime", color = "Branching Path") +
      xlab(xlab) +
      ylab(ylab)
  })
  return(p)
}
