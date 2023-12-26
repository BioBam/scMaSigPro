#' @title Plot trend of the single gene.
#'
#' @description
#' Plot trend of the single gene across the binned pseudotime.
#'
#' @import ggplot2
#' @importFrom RColorConesa getConesaColors
#'
#' @param scmpObj object of class ScMaSigPro.
#' @param feature_id Name of the gene to be plotted.
#' @param xlab X-axis label. (Default is "Pooled Pseudotime")
#' @param ylab Y-axis label. (Default is "Pseudobulk Expression")
#' @param smoothness How smooth the trend should be. Setting to
#' higher values will result in more linear trends. (Default is 0.01)
#' @param logs Whether to log transform counts. (Default is TRUE)
#' @param logType How to log transform the values. Available options 'log',
#' 'log2', 'log10'. (Default is 'log')
#' @param pseudoCount Add a pseudo-count before taking the log. (Default is 1)
#' @param significant Plot gene only if the models are significant based on
#' \code{\link{scMaSigPro::sc.filter()}}. (Default is TRUE)
#' @param summary_mode Compress the expression values for per replicate
#' (if present) per binned pseudotime point. Default is 'median'. Other option
#' 'mean'
#'
#' @return ggplot2 plot object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
plotTrend <- function(scmpObj,
                      feature_id,
                      xlab = "Pooled Pseudotime",
                      ylab = "Pseudobulk Expression",
                      smoothness = 0.01,
                      logs = TRUE,
                      logType = "log",
                      pseudoCount = 1,
                      significant = TRUE,
                      summary_mode = "median") {
  # Invoke Variables
  pb.counts <- "pb.counts"
  pooled.time <- "pooled.time"
  path <- "path"

  # Check summary_mode
  assert_that(any(summary_mode %in% c("median", "mean")),
    msg = paste(
      paste0("'", summary_mode, "'"), "is not a valid option. Please use one of",
      paste(c("median", "mean"), collapse = ", ")
    )
  )

  # Extract edisgn
  alloc.frame <- scmpObj@Design@assignment_matrix %>% as.data.frame()

  # Extract the bulk counts
  bulk.counts <- scmpObj@Dense@assays@data@listData$bulk.counts

  # Check
  assert_that(all(feature_id %in% rownames(bulk.counts)),
    msg = "Feature Id doesn't exist please select another one"
  )

  if (significant) {
    assert_that(any(feature_id %in% unique(unlist(scmpObj@Significant@genes))),
      msg = "Feature Id didn't pass the R2 threshold, please re-run sc.filter, with lower a value or set 'significant' to 'FALSE'"
    )
  }

  # gene_i
  yy <- bulk.counts[rownames(bulk.counts) %in% feature_id, , drop = FALSE]

  # group Vector
  groups.vector <- scmpObj@Design@groups.vector

  # Prepare for Tfit
  rm <- matrix(yy, nrow = 1, ncol = length(yy))
  rownames(rm) <- c("ratio medio")
  colnames(rm) <- rownames(scmpObj@Design@predictor_matrix)

  # Extract the beta
  betas.table <- showCoeff(scmpObj, view = FALSE, return = TRUE)
  betas <- betas.table[rownames(betas.table) %in% feature_id, , drop = FALSE]

  # Set Data
  curve.df <- data.frame(x = 0, y = 0, path = scmpObj@Parameters@path_prefix)
  line.df <- data.frame(x = 0, y = 0, path = scmpObj@Parameters@path_prefix)
  colnames(line.df) <- c("x", "y", scmpObj@Parameters@path_col)
  colnames(line.df) <- c("x", "y", scmpObj@Parameters@path_col)
  curve_data <- NULL
  path.names <- unique(scmpObj@Dense@colData[[scmpObj@Parameters@path_col]])

  # Get x and y
  x <- y <- rep(0, nrow(alloc.frame))

  # Create Point df
  points.df <- data.frame(
    pooled.time = alloc.frame[, scmpObj@Parameters@bin_ptime_col],
    pb.counts = as.vector(t(as.matrix(yy))),
    path = scmpObj@Dense@colData[[scmpObj@Parameters@path_col]]
  )
  # View(yy)
  # View(points.df)
  # stop()
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
    time <- alloc.frame[alloc.frame[[i]] == 1, scmpObj@Parameters@bin_ptime_col]

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

  # Generate line.df
  line.df <- points.df

  # Apply Summary Operation
  if (summary_mode == "mean") {
    line.df <- line.df %>%
      group_by(pooled.time, path) %>%
      summarize(pb.counts = mean(pb.counts), .groups = "drop")
  } else if (summary_mode == "median") {
    line.df <- line.df %>%
      group_by(pooled.time, path) %>%
      summarize(pb.counts = median(pb.counts), .groups = "drop")
  }

  # Plot
  p <- ggplot() +
    geom_point(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), fill = "#102C57", alpha = 0.5, size = 2, stroke = 1, shape = 21) +
    geom_line(data = line.df, aes(x = pooled.time, y = .data$pb.counts, color = path), linetype = "solid", linewidth = 1, alpha = 0.7) +
    geom_line(data = curve.df, aes(x = x, y = y, color = path), linetype = "dashed", linewidth = 1, alpha = 0.7) +
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
