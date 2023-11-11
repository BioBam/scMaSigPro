#' Plot Groups Function
#'
#' This function generates plots based on various parameters. It calculates the summary mode, colors, and other visual attributes to create a plot.
#'
#' @param scmpObj object of class scmpObj
#' @param groupBy A vector of maximum length 6 or a single feature of ID
#' @param xlab X-axis label. Default is "Pooled Pseudotime".
#' @param ylab Y-axis label. Default is "Pseudobulk Expression".
#' @param smoothness abc
#' @param logs Whether to plot log of counts
#' @param logType Log type required
#'
#' @import ggplot2
#' @importFrom RColorConesa getConesaColors
#' @return Generates a plot.
#' @export

sc.PlotProfiles <-
  function(scmpObj,
           xlab = "Pooled Pseudotime",
           ylab = "Pseudobulk Expression",
           smoothness = 0.01,
           logs = TRUE,
           logType = "log10",
           groupBy = "coeff") {
    # Invoke Variables
    pb.counts <- "pb.counts"
    pooled.time <- "pooled.time"
    path <- "path"

    # Extract edisgn
    edesign.frame <- scmpObj@edesign@edesign %>% as.data.frame()

    # Extract the bulk counts
    bulk.counts <- scmpObj@compress.sce@assays@data@listData$bulk.counts

    if (groupBy == "coeff") {
      cluster_vector <- scmpObj@sig.genes@feature.clusters[["sigCoeff"]]
      yy_mat <- bulk.counts[rownames(bulk.counts) %in% names(
        cluster_vector
      ), , drop = FALSE]
    } else if (groupBy == "feature") {
      cluster_vector <- scmpObj@sig.genes@feature.clusters[["sigCounts"]]
      yy_mat <- bulk.counts[rownames(bulk.counts) %in% names(
        cluster_vector
      ), , drop = FALSE]
    }

    # Extract the bulk counts
    edesign <- edesign.frame

    # group Vector
    groups_vector <- scmpObj@scPVector@groups.vector

    # Generate Trends
    trend.data.list <- lapply(c(1:nrow(yy_mat)), function(i, groups.vector = groups_vector, cluster.vector = cluster_vector,
                                                          group_by = groupBy) {
      # Generate vector
      yy <- unlist(yy_mat[i, , drop = TRUE])

      # print(yy)

      # Get featureid
      feature_id <- rownames(yy_mat[i, , drop = FALSE])

      # print(feature_id)

      # Prepare for Tfit
      rm <- matrix(yy, nrow = 1, ncol = length(yy))
      rownames(rm) <- c("ratio medio")
      colnames(rm) <- rownames(scmpObj@edesign@dis)

      # print(rm)

      # Extract the beta
      betas.table <- showCoeff(scmpObj, view = FALSE, return = TRUE)
      betas <- betas.table[rownames(betas.table) %in% feature_id, , drop = FALSE]

      # print(betas)

      # Set Data
      curve.df <- data.frame(x = 0, y = 0, path = scmpObj@addParams@path_prefix)
      line.df <- data.frame(x = 0, y = 0, path = scmpObj@addParams@path_prefix)
      colnames(line.df) <- c("x", "y", scmpObj@addParams@path_colname)
      colnames(line.df) <- c("x", "y", scmpObj@addParams@path_colname)
      curve_data <- NULL
      path.names <- unique(scmpObj@compress.sce@colData[[scmpObj@addParams@path_colname]])

      # Get x and y
      x <- y <- rep(0, nrow(edesign.frame))

      # Create Point df
      points.df <- data.frame(
        pooled.time = edesign.frame[, scmpObj@addParams@bin_pseudotime_colname],
        pb.counts = as.vector(yy),
        path = scmpObj@compress.sce@colData[[scmpObj@addParams@path_colname]]
      )

      #    print(points.df)

      for (i in path.names) {
        # Extract Coeff
        a <- reg.coeffs(
          coefficients = betas,
          groups.vector = groups.vector,
          group = i
        )
        a <- c(a, rep(0, (7 - length(a))))

        # Extract the time
        time <- edesign.frame[edesign.frame[[i]] == 1, scmpObj@addParams@bin_pseudotime_colname]

        # Create a data frame with time values
        x <- seq(from = min(time), to = max(time), by = smoothness)

        # Compute the curve values
        y <- a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) +
          a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5)

        # Create tmpvector
        curve_df_tmp <- data.frame(
          x = x, y = y,
          path = i
        )
        curve.df <- rbind(curve.df, curve_df_tmp)
      }

      curve.df <- curve.df[-1, ]

      # Calc limits
      xlim <- c(min(points.df[[pooled.time]], na.rm = TRUE), max(points.df[[pooled.time]], na.rm = TRUE) * 1.3)
      # ylim <- c(min(as.numeric(points.df[[pb.counts]]), na.rm = TRUE), max(as.numeric(points.df[[pb.counts]]), na.rm = TRUE))

      xlim[2] <- max(points.df[[pooled.time]])

      conesa_colors <- getConesaColors()[c(TRUE, FALSE)][c(1:length(unique(points.df[[path]])))]
      names(conesa_colors) <- unique(points.df[[path]])

      # Extract sol
      data.sol <- showSol(scmpObj, view = FALSE, return = TRUE)
      data.sol <- data.sol[feature_id, , drop = FALSE]

      # if log is requestion
      if (logs) {
        if (logType == "log2") {
          points.df$pb.counts <- log2(points.df$pb.counts)
          ylab <- paste0("log2(", ylab, ")")
        } else if (logType == "log") {
          points.df$pb.counts <- log(points.df$pb.counts)
          ylab <- paste0("log(", ylab, ")")
        } else if (logType == "log10") {
          points.df$pb.counts <- log10(points.df$pb.counts)
          ylab <- paste0("log10(", ylab, ")")
        } else {
          stop("'logType' should be one of 'log2', 'log10', 'log'")
        }
      }
      if (group_by == "coeff") {
        plot.df <- curve.df
        plot.df[["feature_id"]] <- feature_id
        plot.df[["cluster_id"]] <- cluster.vector[names(cluster.vector) == feature_id]
      } else if (group_by == "feature") {
        plot.df <- points.df
        colnames(plot.df) <- c("x", "y", "path")
        plot.df[["feature_id"]] <- feature_id
        plot.df[["cluster_id"]] <- cluster.vector[names(cluster.vector) == feature_id]
      }

      return(plot.df)
    })

    # Combine
    cluster.trend.data <- do.call("rbind", trend.data.list)

    # Assuming feature_id and cluster_id are factors
    p <- ggplot(data = cluster.trend.data, aes(x = .data$x, y = log(.data$y), group = interaction(.data$feature_id, .data$path), color = .data$path)) +
      geom_line(aes(linetype = .data$path), size = 0.4) + # Draw lines
      geom_point(size = 1, shape = 21, alpha = 0.5, stroke = 1) + # Draw points
      facet_wrap(~ .data$cluster_id, scales = "free_y") + # Create a panel for each cluster_id
      scale_color_manual(values = colorConesa(length(unique(.data$path)))) + # Custom colors for paths
      theme_classic(base_size = 12) +
      # geom_smooth(data = df, aes(x = x, y = log(y), color = path), method = "loess", se = FALSE, size = 0.7)+
      theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, angle = 0),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis text if necessary
      ) +
      labs(title = "Gene Expression over Pseudotime", color = "Path") +
      guides(color = guide_legend(title = "Path")) # Adjust legend for paths
    return(p)
  }
