#' Plot Groups Function
#'
#' This function generates plots based on various parameters. It calculates the summary mode, colors, and other visual attributes to create a plot.
#'
#' @param scmpObj object of class scmpObj
#' @param feature_id A vector of maximum length 6 or a single feature of ID
#' @param xlab X-axis label. Default is "Pooled Pseudotime".
#' @param ylab Y-axis label. Default is "Pseudobulk Expression".
#' @param smoothness abc
#' @param logs Whether to plot log of counts
#' @param logType Log type required
#' @param pseudoCount add pseucount for log
#' @param significant description
#'
#' @import ggplot2
#' @importFrom RColorConesa getConesaColors
#' @return Generates a plot.
#' @export

plotTrend <-
  function(scmpObj,
           feature_id,
           xlab = "Pooled Pseudotime",
           ylab = "Pseudobulk Expression",
           smoothness = 0.01,
           logs = TRUE,
           logType = "log",
           pseudoCount = 1,
           significant = TRUE) {
    # Invoke Variables
    pb.counts <- "pb.counts"
    pooled.time <- "pooled.time"
    path <- "path"

    # Extract edisgn
    alloc.frame <- scmpObj@design@alloc %>% as.data.frame()

    # Extract the bulk counts
    bulk.counts <- scmpObj@dense@assays@data@listData$bulk.counts

    # Check
    assert_that(all(feature_id %in% rownames(bulk.counts)),
      msg = "Feature Id doesn't exist please select another one"
    )

    if (significant) {
      assert_that(any(feature_id %in% unique(unlist(scmpObj@sig.genes@sig.genes))),
        msg = "Feature Id didn't pass the R2 threshold, please re-run sc.get.sigenes, with lower a value or set 'significant' to 'FALSE'"
      )
    }

    # gene_i
    yy <- bulk.counts[rownames(bulk.counts) %in% feature_id, , drop = FALSE]

    # group Vector
    groups.vector <- scmpObj@design@groups.vector

    # Prepare for Tfit
    rm <- matrix(yy, nrow = 1, ncol = length(yy))
    rownames(rm) <- c("ratio medio")
    colnames(rm) <- rownames(scmpObj@design@predictor)

    # Extract the beta
    betas.table <- showCoeff(scmpObj, view = FALSE, return = TRUE)
    betas <- betas.table[rownames(betas.table) %in% feature_id, , drop = FALSE]

    # Set Data
    curve.df <- data.frame(x = 0, y = 0, path = scmpObj@param@path_prefix)
    line.df <- data.frame(x = 0, y = 0, path = scmpObj@param@path_prefix)
    colnames(line.df) <- c("x", "y", scmpObj@param@path_colname)
    colnames(line.df) <- c("x", "y", scmpObj@param@path_colname)
    curve_data <- NULL
    path.names <- unique(scmpObj@dense@colData[[scmpObj@param@path_colname]])

    # Get x and y
    x <- y <- rep(0, nrow(alloc.frame))

    # Create Point df
    points.df <- data.frame(
      pooled.time = alloc.frame[, scmpObj@param@bin_pseudotime_colname],
      pb.counts = as.vector(t(yy)),
      path = scmpObj@dense@colData[[scmpObj@param@path_colname]]
    )
    #View(yy)
    #View(points.df)
    #stop()
    for (i in path.names) {
      # Extract Coeff
      a <- reg.coeffs(
        coefficients = betas,
        groups.vector = groups.vector,
        group = i
      )
      a <- c(a, rep(0, (7 - length(a))))
      a[is.na(a)] <- 0

      # Extract the time
      time <- alloc.frame[alloc.frame[[i]] == 1, scmpObj@param@bin_pseudotime_colname]

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
        points.df$pb.counts <- log2(points.df$pb.counts + pseudoCount)
        ylab <- paste0("log2(", ylab, ")")
      } else if (logType == "log") {
        points.df$pb.counts <- log(points.df$pb.counts + pseudoCount)
        ylab <- paste0("log(", ylab, ")")
      } else if (logType == "log10") {
        points.df$pb.counts <- log10(points.df$pb.counts + pseudoCount)
        ylab <- paste0("log10(", ylab, ")")
      } else {
        stop("'logType' should be one of 'log2', 'log10', 'log'")
      }
    }

    # Convert Negative values to 0
    #curve.df[curve.df$y < 0, "y"] <- 0
    # points.df[points.df$pb.counts < 0, "pb.counts"] <- 0
    
    # Generate line.df
    line.df <- points.df
    
    line.df <- line.df %>%
        group_by(pooled.time, path) %>%
        summarize(mean_pb_counts = mean(pb.counts),  .groups = "drop")
    
    p <- ggplot() +
      geom_point(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), fill = "#102C57", alpha = 0.5, size = 2, stroke = 1, shape = 21) +
      geom_line(data = line.df, aes(x = pooled.time, y = mean_pb_counts, color = path), linetype = "dotted", linewidth = 1) +
      geom_line(data = curve.df, aes(x = x, y = y, color = path), linetype = "solid", linewidth = 1.5) +
      ggtitle(
        paste("Feature Id:", feature_id),
        subtitle = paste("R2:", round(data.sol[, 2], 3), "| p-Value:", round(data.sol[, 1], 3))
      ) +
      xlab(xlab) +
      ylab(ylab) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = seq(min(xlim), max(xlim), by = round(log10(length(points.df[[pooled.time]]))))) +
      labs(color = "Paths") +
      # coord_cartesian(xlim = xlim, ylim = ylim) +
      scale_color_manual(values = conesa_colors)
    #
    return(p)
  }
