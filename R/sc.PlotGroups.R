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
#' @param ... Additional plotting parameters.
#'
#' @return Generates a plot.
#' @export

sc.PlotGroups <-
  function(scmpObj, feature_id, edesign = NULL, time = edesign[, 1], groups = edesign[, c(3:ncol(edesign))], repvect = edesign[, 2],
           dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05, nvar.correction = FALSE,
           summary.mode = "median", groups.vector = NULL, main = NULL, sub = NULL,
           xlab = "Pooled Pseudotime", ylab = "Pseudobulk Expression", item = NULL, ylim = NULL, pch = 21,
           col = NULL, legend = TRUE, cex.legend = 1, show.umap = T,
           lty.legend = NULL, ...) {
      
      
      data <- scmpObj@scTFit@dat
      assert_that(all(feature_id %in% rownames(data)),
                    msg = "Feature Id doesn't exist please select another one")
      data <- data[rownames(data) %in% feature_id, , drop = F]
      
      sol <- showSol(scmpObj, view = F, return = T)
      data.sol <- sol[rownames(sol) %in% feature_id, , drop = F]
      # Get the R2 and Pvalues
      
      
      #View(data)
      
    # Check if data is a vector
    if (!is.vector(data)) {
      if (summary.mode == "representative") {
        distances <- apply(as.matrix(dist(data,
          diag = TRUE,
          upper = TRUE
        )), 1, sum)
        representative <- names(distances)[distances == min(distances)]
        yy <- as.numeric(data[rownames(data) == representative, ])
        sub <- paste("Representative:", representative)
      } else if (summary.mode == "median") {
        yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
        if (is.null(sub)) {
          sub <- paste("Median profile of", nrow(data), item, sep = " ")
        }
      } else {
        stop("not valid summary.mode")
      }
      if (dim(data)[1] == 1) {
        main <- rownames(data)
        sub <- NULL
      }
    } else if (length(data) != 0) {
      yy <- as.numeric(data)
      sub <- rownames(data)
    } else {
      stop("empty data")
    }
    if (is.null(ncol(groups))) {
      ncol <- 1
      legend <- FALSE
      codeg <- "group"
    } else {
      ncol <- ncol(groups)
      codeg <- as.character(colnames(groups))
    }

    reps <- i.rank(repvect)
    y <- vector(mode = "numeric", length = length(unique(reps)))
    x <- vector(mode = "numeric", length = length(unique(reps)))
    g <- matrix(nrow = length(unique(reps)), ncol = ncol)
    for (k in 1:length(y)) {
      y[k] <- mean(yy[reps == k], na.rm = TRUE)
      x[k] <- mean(time[reps == k])
      for (j in 1:ncol) {
        g[k, j] <- mean(groups[reps == k, j])
      }
    }
    if (is.null(ylim)) {
      ylim <- c(min(as.numeric(yy), na.rm = TRUE), max(as.numeric(yy),
        na.rm = TRUE
      ))
    }
    abcissa <- x
    xlim <- c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE) *
      1.3)
    if (is.null(col)) {
      color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
      color2 <- groups
      for (j in 1:ncol) {
        color2[, j] <- color2[, j] * j
      }
      color2 <- as.vector(apply(color2, 1, sum) + 1)
    } else {
      color1 <- col
      color2 <- groups
      for (j in 1:ncol) {
        color2[, j] <- color2[, j] * col[j]
      }
      color2 <- as.vector(apply(color2, 1, sum))
    }

    paths <- setNames(colnames(groups), unique(color2))
    paths <- paths[as.character(color2)]
    names(paths) <- NULL

    points.df <- data.frame(
      pooled.time = time,
      pb.counts = yy,
      path = paths
    )
    rm <- matrix(yy, nrow = 1, ncol = length(yy))
    rownames(rm) <- c("ratio medio")
    colnames(rm) <- rownames(dis)
    fit.y <- T.fit(rm,
      design = dis, step.method = step.method,
      min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction
    )
    betas <- fit.y$coefficients

    curve.df <- data.frame(x = 0, y = 0, path = "path")
    line.df <- data.frame(x = 0, y = 0, path = "path")
    curve_data <- NULL
    path.names <- colnames(groups)

    for (i in 1:ncol(groups)) {
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
      lines(lxo, lyo, col = color1[i], ...)

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
              subtitle = paste("R2:", round(data.sol[,2], 3), "| p-Value:", round(data.sol[,1], 3))) +
      xlab(xlab) +
      ylab(ylab) +
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom",
            panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
            panel.grid.minor = element_blank()) +
        scale_x_continuous(breaks = seq(min(xlim), max(xlim), by = round(log10(length(time)))))+
      labs(color = "Paths") +
        coord_cartesian(xlim = xlim, ylim = ylim)+
    
      scale_color_manual(values = conesa_colors)


    print(p)

  }
