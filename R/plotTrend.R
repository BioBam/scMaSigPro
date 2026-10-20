#' @title Plot trend of the single gene.
#'
#' @description
#' Plot trend of the single gene across the binned pseudotime.
#'
#' @importFrom ggplot2 aes element_blank element_line geom_line geom_point
#'   ggplot ggtitle labs scale_color_manual scale_x_continuous theme
#'   theme_classic xlab ylab
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
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
#' \code{scMaSigPro::sc.filter()}. (Default is TRUE)
#' @param summary_mode Compress the expression values per replicate (if present)
#'  per binned pseudotime point. Default is 'median'. Other option 'mean'
#' @param curves Whether to plot the fitted curves. (Default is TRUE)
#' @param lines Whether to plot the lines. (Default is FALSE)
#' @param points Whether to plot the points. (Default is TRUE)
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
                      summary_mode = "median",
                      curves = TRUE,
                      lines = FALSE,
                      points = TRUE) {
    # Invoke Variables
    pb.counts <- "pb.counts"
    pooled.time <- "pooled.time"
    path <- "path"
    x <- "x"
    y <- "y"

    # Check Object Validity
    .check_scmp(scmpObj)

    # Check summary_mode
    summary_mode <- .check_choice(summary_mode, c("median", "mean"),
        msg = paste(
            paste0("'", summary_mode, "'"),
            "is not a valid option. Please use one of",
            paste(c("median", "mean"), collapse = ", ")
        )
    )

    # Check Assertion
    .check(curves || lines || points,
        msg = "At least one of 'curves', 'lines', or 'points' must be TRUE."
    )
    if (logs) {
        logType <- .check_choice(logType, c("log", "log2", "log10"),
            msg = "'logType' should be one of 'log2', 'log10', 'log'"
        )
    }

    # The data of the layers
    trend <- .trendData(scmpObj, feature_id,
        smoothness = smoothness, logs = logs, logType = logType,
        pseudoCount = pseudoCount, significant = significant,
        summary_mode = summary_mode
    )
    points.df <- trend$points
    line.df <- trend$lines
    curve.df <- trend$curves
    if (logs) {
        ylab <- paste0(logType, "(", ylab, ")")
    }

    # Calc limits
    xlim <- c(min(points.df[[pooled.time]], na.rm = TRUE), max(points.df[[pooled.time]], na.rm = TRUE) * 1.3)

    xlim[2] <- max(points.df[[pooled.time]])

    scmp_pal <- scmp_colors(n = length(trend$path.names))
    names(scmp_pal) <- unique(points.df[[path]])

    # Extract sol
    data.sol <- showSol(scmpObj)
    data.sol <- data.sol[feature_id, , drop = FALSE]

    # Plot
    layer_names <- c(NULL)
    p <- ggplot() +
        ggtitle(
            paste("Feature Id:", feature_id),
            subtitle = paste("R2:", round(data.sol[, 2], 3), "| p-Value:", round(data.sol[, 1], 3))
        ) +
        xlab(xlab) +
        ylab(ylab)
    names(p$layers) <- layer_names

    if (points) {
        p <- p + geom_point(data = points.df, aes(x = pooled.time, y = pb.counts, color = path), fill = "#102C57", alpha = 0.4, size = 1.5, stroke = 1, shape = 21)
        layer_names <- c(layer_names, "points")
        names(p$layers) <- layer_names
    }
    if (lines) {
        p <- p + geom_line(data = line.df, aes(x = pooled.time, y = pb.counts, color = path), linetype = "dashed", linewidth = 0.6, alpha = 0.7)
        layer_names <- c(layer_names, "lines")
        names(p$layers) <- layer_names
    }
    if (curves) {
        p <- p + geom_line(data = curve.df, aes(x = x, y = y, color = path), linetype = "solid", linewidth = 0.7, alpha = 0.8)
        layer_names <- c(layer_names, "curves")
        names(p$layers) <- layer_names
    }


    p <- p + theme_classic(base_size = 12) +
        theme(
            legend.position = "bottom",
            panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
            panel.grid.minor = element_blank()
        ) +
        scale_x_continuous(breaks = seq(min(xlim), max(xlim), by = round(log10(length(points.df[[pooled.time]]))))) +
        labs(color = "Paths") +
        # coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_color_manual(values = scmp_pal)
    return(p)
}

# The data of the layers of plotTrend() for one gene: the points
# (pooled.time, pb.counts, path), the lines (the points summarised per
# binned time and path, or the points themselves with an offset) and the
# curve of the fitted polynomial of each path (x, y, path), with the path
# names. plotTrendCluster() collects the same frames for the genes of each
# cluster.
.trendData <- function(scmpObj, feature_id, smoothness, logs, logType,
                       pseudoCount, significant, summary_mode) {
    # Offset
    offset_vector <- scmpObj@Design@offset

    # Extract edisgn
    alloc.frame <- as.data.frame(scmpObj@Design@assignment_matrix)

    # Extract the bulk counts
    .check_dense_counts(scmpObj)
    bulk.counts <- assay(scmpObj@Dense, "bulk.counts", withDimnames = FALSE)

    # Check
    .check(all(feature_id %in% rownames(bulk.counts)),
        msg = paste0(
            "'", feature_id, "' doesn't exist please select another one."
        )
    )

    if (significant) {
        .check(
            any(feature_id %in% unique(unlist(scmpObj@Significant@genes))),
            msg = paste(
                "Feature Id didn't pass the R2 threshold, please re-run",
                "sc.filter() with a lower 'rsq' value or set 'significant' to",
                "'FALSE'"
            )
        )
    }

    # gene_i
    yy <- bulk.counts[rownames(bulk.counts) %in% feature_id, , drop = FALSE]

    # group Vector
    groups.vector <- scmpObj@Design@groups.vector

    # Extract the beta
    betas.table <- showCoeff(scmpObj)
    betas <- betas.table[rownames(betas.table) %in% feature_id, , drop = FALSE]

    # Set Data
    curve.df <- data.frame(x = 0, y = 0, path = scmpObj@Parameters@path_prefix)
    path.names <- unique(
        colData(scmpObj@Dense)[[scmpObj@Parameters@path_col]]
    )

    # Create Point df
    points.df <- data.frame(
        pooled.time = alloc.frame[, scmpObj@Parameters@bin_ptime_col],
        pb.counts = as.vector(t(as.matrix(yy))),
        path = colData(scmpObj@Dense)[[scmpObj@Parameters@path_col]]
    )
    for (i in path.names) {
        # Extract Coeff
        a <- .path_coeffs(
            coefficients = betas,
            groups.vector = groups.vector,
            group = i
        )
        a[is.na(a)] <- 0

        # Extract the time
        time <- alloc.frame[alloc.frame[[i]] == 1, scmpObj@Parameters@bin_ptime_col]

        # Create a data frame with time values
        x <- seq(from = min(time), to = max(time), by = smoothness)

        # The polynomial of the path, one power of x per coefficient, added
        # in order: a[1] + a[2] * x + a[3] * x^2 + ...
        y <- a[1] + a[2] * x
        for (k in seq_along(a)[-c(1L, 2L)]) {
            y <- y + a[k] * (x^(k - 1))
        }

        # Create tmpvector
        curve_df_tmp <- data.frame(
            x = x, y = y,
            path = i
        )
        curve.df <- rbind(curve.df, curve_df_tmp)
    }

    curve.df <- curve.df[-1, ]

    # Correct by offset
    if (sum(offset_vector) != 0) {
        points.df["pb.counts"] <- points.df["pb.counts"] / exp(offset_vector)
    }

    # if log is requestion
    suppressWarnings(
        expr = {
            if (logs) {
                if (logType == "log2") {
                    points.df$pb.counts <- log2(points.df$pb.counts + pseudoCount)
                } else if (logType == "log") {
                    points.df$pb.counts <- log(points.df$pb.counts + pseudoCount)
                } else if (logType == "log10") {
                    points.df$pb.counts <- log10(points.df$pb.counts + pseudoCount)
                }
            }
        }
    )

    # Generate line.df
    line.df <- points.df

    # Apply Summary Operation: one value per binned time and path
    if (summary_mode == "mean") {
        line.df <- .summarise_by(
            line.df, c("pooled.time", "path"), "pb.counts", mean
        )
    } else if (summary_mode == "median") {
        line.df <- .summarise_by(
            line.df, c("pooled.time", "path"), "pb.counts", median
        )
    }

    if (sum(offset_vector) != 0) {
        line.df <- points.df
    }

    list(
        points = points.df, lines = line.df, curves = curve.df,
        path.names = path.names
    )
}
