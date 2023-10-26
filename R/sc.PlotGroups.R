#' Plot Groups Function
#'
#' This function generates plots based on various parameters. It calculates the summary mode, colors, and other visual attributes to create a plot.
#'
#' @param scmpObj object of class scmpObj
#' @param feature_id A vector of maximum length 6 or a single feature of ID
#' @param edesign Experimental design matrix.
#' @param time Vector indicating time points.
#' @param groups A matrix representing the groups.
#' @param repvect Repetition vector.
#' @param dis Distance matrix.
#' @param step.method Method for backward elimination. Default is "backward".
#' @param min.obs Minimum number of observations. Default is 2.
#' @param alfa Significance level for statistical tests. Default is 0.05.
#' @param nvar.correction Logical, should variance correction be applied? Default is FALSE.
#' @param summary.mode Mode of summary. Either "representative" or "median". Default is "median".
#' @param groups.vector Additional group vectors.
#' @param sub Subtitle of the plot.
#' @param xlab X-axis label. Default is "Pooled Pseudotime".
#' @param ylab Y-axis label. Default is "Pseudobulk Expression".
#' @param item Optional items to be included.
#' @param ylim Y-axis limits.
#' @param pch Plotting symbol to be used.
#' @param col Colors for plot elements.
#' @param legend Logical, should legend be displayed? Default is TRUE.
#' @param cex.legend Size of legend text. Default is 1.
#' @param lty.legend Line type for the legend.
#' @param plot.bins Additionally Plot bins
#' @param ... Additional plotting parameters.
#'
#' @importFrom maSigPro i.rank T.fit reg.coeffs
#' @importFrom ggplot2 ggplot
#' @return Generates a plot.
#' @export

sc.PlotGroups <-
  function(scmpObj, feature_id, edesign = NULL, time = edesign[, 1], groups = edesign[, c(3:ncol(edesign))], repvect = edesign[, 2],
           dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05, nvar.correction = FALSE,
           summary.mode = "median", groups.vector = NULL, main = NULL, sub = NULL,
           xlab = "Pooled Pseudotime", ylab = "Pseudobulk Expression", item = NULL, ylim = NULL, pch = 21,
           col = NULL, legend = TRUE, cex.legend = 1, show.umap = T,
           lty.legend = NULL, plot.bins = FALSE) {
      
      # Extract the bulk counts
      bulk.counts = scmpObj@compress.sce@assays@data@listData$bulk.counts
      
      # Check
      assert_that(all(feature_id %in% rownames(bulk.counts)),
                  msg = "Feature Id doesn't exist please select another one"
      )
      # gene_i
      yy <- bulk.counts[rownames(bulk.counts) %in% feature_id, , drop = F]
      
      # Extract the bulk counts
      edesign = scmpObj@edesign@edesign
      
      # group Vector
      groups.vector = scmpObj@scPVector@groups.vector
      
      # Prepare for Tfit
      rm <- matrix(yy, nrow = 1, ncol = length(yy))
      rownames(rm) <- c("ratio medio")
      colnames(rm) <- rownames(scmpObj@edesign@dis)
      
      # Extract the beta
      betas.table <- showCoeff(scmpObj, view = F, return = T)
      betas <- betas.table[feature_id, , drop = F]
      
      # Set Data
      curve.df <- data.frame(x = 0, y = 0, path = scmpObj@addParams@path_prefix)
      line.df <- data.frame(x = 0, y = 0, path = scmpObj@addParams@path_prefix)
      colnames(line.df) <- c("x", "y", scmpObj@addParams@path_colname)
      colnames(line.df) <- c("x", "y", scmpObj@addParams@path_colname)
      curve_data <- NULL
      path.names <- unique(scmpObj@compress.sce@colData[[scmp@addParams@path_colname]])
      
      # Get x and y
      x <- y <- rep(0, nrow(scmpObj@edesign@edesign))
      
      PlotGroups(data = yy,
                 edesign = scmp@edesign@edesign,
                 show.lines = T,
                 show.fit = T,
                 dis = scmp@edesign@dis,
                 groups.vector = scmp@scPVector@groups.vector,
                 summary.mode = "median"
                 
      )
      
      # Create Point df
      points.df <- data.frame(
          pooled.time = scmpObj@edesign@edesign[, scmpObj@addParams@bin_pseudotime_colname],
          pb.counts = as.vector(yy),
          path = scmpObj@compress.sce@colData[[scmpObj@addParams@path_colname]]
      )
      

      
      for (i in path.names){
          
          # Extract Coeff
          a <- reg.coeffs(
              coefficients = betas, 
              groups.vector = groups.vector,
              group = i
          )
          a <- c(a, rep(0, (7 - length(a))))
          
          # Extract the time
          time <- scmpObj@edesign@edesign[scmpObj@edesign@edesign[[i]] == 1, scmpObj@addParams@bin_pseudotime_colname]
          
          # Create a data frame with time values
          x <- seq(from = min(time), to = max(time), by = 0.001)
          
          # Compute the curve values
          y <- a[1] + a[2]*x + a[3]*(x^2) + a[4]*(x^3) +
              a[5]*(x^4) + a[6]*(x^5) + a[7]*(x^5)
          
          # Create tmpvector
          curve_df_tmp <- data.frame(
              x = x, y = y,
              path = i
          )
          curve.df <- rbind(curve.df, curve_df_tmp)
      }
      
      curve.df <- curve.df[-1, ]
      
      View(curve.df)
      
      
      
      PlotGroups(data = yy,
                 edesign = scmp@edesign@edesign,
                 show.lines = T,
                 show.fit = T,
                 dis = scmp@edesign@dis,
                 groups.vector = scmp@scPVector@groups.vector,
                 summary.mode = "median"
                 
      )
      
      # Calc limits
      xlim <- c(min(points.df$pooled.time, na.rm = TRUE), max(points.df$pooled.time, na.rm = TRUE) * 1.3)
      ylim <- c(min(as.numeric(points.df$pb.counts), na.rm = TRUE), max(as.numeric(points.df$pb.counts), na.rm = TRUE))
      
      xlim[2] <- max(points.df$pooled.time)
      
      conesa_colors <- getConesaColors()[c(T, F)][c(1:length(unique(points.df$path)))]
      names(conesa_colors) <- unique(points.df$path)
      
      print(xlim)
      p <- ggplot() +
          geom_point(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), fill = "#102C57", alpha = 0.5, size = 2, stroke = 1, shape = 21) +
          geom_line(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), linetype = "dotted", linewidth = 1) +
          geom_line(data = curve.df, aes(x = x, y = y, color = path), linetype = "solid", linewidth = 1.5) +
          ggtitle(paste("Feature Id:", feature_id),
               #   subtitle = paste("R2:", round(data.sol[, 2], 3), "| p-Value:", round(data.sol[, 1], 3))
          ) +
          # xlab(xlab) +
          # ylab(ylab) +
          theme_classic(base_size = 12) +
          theme(
              legend.position = "bottom",
              panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
              panel.grid.minor = element_blank()
          ) +
          scale_x_continuous(breaks = seq(min(xlim), max(xlim), by = round(log10(length(points.df$pooled.time))))) +
          labs(color = "Paths") +
          coord_cartesian(xlim = xlim, ylim = ylim) +
          scale_color_manual(values = conesa_colors)
          # 
      print(p)
      
      stop()
      
      points.df$path <- as.factor(points.df$path)
      curve.df$path <- as.factor(curve.df$path)
      line.df$path <- as.factor(line.df$path)
      
      conesa_colors <- getConesaColors()[c(T, F)][c(1:length(unique(paths)))]
      names(conesa_colors) <- unique(paths)
      
      
      p <- ggplot() +
          geom_point(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), fill = "#102C57", alpha = 0.5, size = 2, stroke = 1, shape = 21) +
          geom_line(data = curve.df, aes(x = x, y = y, color = path), linetype = "solid", linewidth = 1.5) +
          geom_line(data = line.df, aes(x = x, y = y, color = path), linetype = "dotted", linewidth = 1) +
          ggtitle(paste("Feature Id:", feature_id),
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
          scale_x_continuous(breaks = seq(min(xlim), max(xlim), by = round(log10(length(time))))) +
          labs(color = "Paths") +
          coord_cartesian(xlim = xlim, ylim = ylim) +
          scale_color_manual(values = conesa_colors)
      
      
      stop()
      
      # Extract the curve
      for (i in groups) {
          group <- g[, i]
          if (!is.null(betas)) {
              li <- c(2:6)
              a <- reg.coeffs(
                  coefficients = betas, groups.vector = groups.vector,
                  group = colnames(groups)[i]
              )
              
              a <- c(a, rep(0, (7 - length(a))))
              
              # Create a null device to suppress plotting
              invisible({
                  curve_data <- curve(
                      a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) +
                          a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5),
                      from = min(time), to = max(time),
                      col = color1[i], add = F, lty = li[i]
                  )
              })
          }
          
          curve_df_tmp <- data.frame(
              x = curve_data$x, y = curve_data$y,
              path = path.names[i]
          )
          
          curve.df <- rbind(curve.df, curve_df_tmp)
          
          lx <- abcissa[group != 0]
          ly <- y[group != 0]
          ord <- order(lx)
          lxo <- lx[ord]
          lyo <- ly[ord]
          # lines(lxo, lyo, col = color1[i], ...)
          
          line_df_tmp <- data.frame(
              x = lxo, y = lyo,
              path = path.names[i]
          )
          
          line.df <- rbind(line.df, line_df_tmp)
      }
      
      # Factor Convert
      curve.df <- curve.df[-1, ]
      line.df <- line.df[-1, ]
      


    points.df$path <- as.factor(points.df$path)
    curve.df$path <- as.factor(curve.df$path)
    line.df$path <- as.factor(line.df$path)

    conesa_colors <- getConesaColors()[c(T, F)][c(1:length(unique(paths)))]
    names(conesa_colors) <- unique(paths)

    p <- ggplot() +
      geom_point(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), fill = "#102C57", alpha = 0.5, size = 2, stroke = 1, shape = 21) +
      geom_line(data = curve.df, aes(x = x, y = y, color = path), linetype = "solid", linewidth = 1.5) +
      geom_line(data = line.df, aes(x = x, y = y, color = path), linetype = "dotted", linewidth = 1) +
      ggtitle(paste("Feature Id:", feature_id),
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
      scale_x_continuous(breaks = seq(min(xlim), max(xlim), by = round(log10(length(time))))) +
      labs(color = "Paths") +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      scale_color_manual(values = conesa_colors)


    # if(plot.bins){
    #
    #     raw.data <- scmp.obj@sce@assays@data@listData$counts
    #     raw.data.gene <- raw.data[rownames(raw.data) %in% feature_id,,drop = F]
    #     raw.data.gene <- as.data.frame(t(raw.data.gene))
    #     colnames(raw.data.gene) <- "raw_count"
    #     raw.data.gene$cell <- rownames(raw.data.gene)
    #
    #
    #     # Get bins
    #     compressed.df <- as.data.frame(SingleCellExperiment::colData(scmp.obj@compress.sce))
    #
    #     # Select the columns
    #     compressed.df <- compressed.df %>%
    #         separate_rows(cluster.members, sep = "\\|") %>%
    #         as.data.frame()
    #
    #     expand.df <- as.data.frame(SingleCellExperiment::colData(scmp.obj@sce))
    #
    #
    #     # Merge
    #     plt.df <- merge(compressed.df, raw.data.gene, by.x = "cluster.members", by.y = "cell")
    #     plt.df$sd <- sd(raw.data.gene$raw_count)
    #     plt.df$mean <- mean(raw.data.gene$raw_count)
    #     plt.df$sum <- sum(raw.data.gene$raw_count)
    #     plt.df <- merge(plt.df, expand.df, by.x = "cluster.members", by.y = "Cell")
    #
    #
    #     agg_data <- df %>%
    #         group_by(binnedTime, path, bin, bin.size) %>%
    #         summarise(
    #             sum_raw_count = mean(raw_count),
    #             sd_raw_count = sd(raw_count)
    #         )
    #
    #     # Make sure to remove NA (if any) before plotting
    #     agg_data <- na.omit(agg_data)
    #
    #     # Plotting
    #     ggplot(agg_data, aes(x = factor(binnedTime), y = sum_raw_count, color = path)) +
    #         geom_point(size = 3) +
    #         geom_errorbar(aes(ymin = sum_raw_count - sd_raw_count, ymax = sum_raw_count + sd_raw_count), width = 0.2) +
    #         facet_wrap(~path, scales = "free") +
    #         labs(
    #             title = "Aggregated Raw Counts Along the Step",
    #             x = "Step",
    #             y = "Aggregated Raw Counts"
    #         ) +
    #         theme_minimal() +
    #         theme(legend.position = "bottom")
    #
    #
    #     ggplot(data = df,
    #            aes(axis1 = Step, axis2 = binnedTime, y = raw_count)) +
    #         geom_alluvium(aes(fill = binnedTime)) +
    #         geom_stratum() +
    #         facet_wrap(~path, scales = "free") +
    #         geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    #         labs(title = "Sankey Diagram of Step and Binned Time") +
    #         theme_minimal()
    #
    #
    # }
    print(p)
  }
