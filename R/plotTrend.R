#' @title Plot the trend of a gene
#'
#' @description
#' `plotTrend()` plots the pseudo-bulk counts of a gene along the binned
#' pseudotime of each path, with the curves of its selected model
#' ([sc.t.fit()]).
#'
#' @param xlab The label of the x axis.
#' @param ylab The label of the y axis; with logged points it becomes, for
#'   example, `log(ylab)`.
#' @param smoothness The step of the pseudotime grid of the curves: the
#'   smaller, the smoother.
#' @param logs Whether to log the points: the counts of the bins, divided,
#'   for log-link fits, by the size factors stored in the design. The fitted
#'   curve is the polynomial of the selected model on the scale of the link:
#'   the log of the mean under the log link, the mean itself under the
#'   identity link. The default, `NULL`, logs the points when the fits used
#'   the log link and draws them as they are otherwise, so that points and
#'   curve are on the same scale.
#' @param logType The log of the points: `"log"`, `"log2"` or `"log10"`.
#' @param pseudoCount The number added to the points before the log.
#' @param significant Whether the gene must be in a gene list of
#'   [sc.filter()].
#' @param summary_mode How the points of each binned pseudotime and path are
#'   summarised for the lines: `"median"` or `"mean"`.
#' @param curves Whether to draw the fitted curves.
#' @param lines Whether to draw lines through the points.
#' @param points Whether to draw the points.
#'
#' @importFrom ggplot2 aes element_blank element_line geom_line geom_point
#'   ggplot ggtitle labs scale_color_manual scale_x_continuous theme
#'   theme_classic xlab ylab
#'
#' @return A ggplot object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @family plots
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' plotTrend(scmp.ob, "Gene10")
#'
#' @export
plotTrend <- function(scmpObj,
                      feature_id,
                      xlab = "Pooled Pseudotime",
                      ylab = "Pseudobulk Expression",
                      smoothness = 0.01,
                      logs = NULL,
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
    logs <- .trend_logs(scmpObj, logs)

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

# The default of `logs` in plotTrend() and plotTrendCluster():
# the curve is the polynomial on the scale of the link, so the points are
# logged when the fits used the log link and left as they are otherwise.
# The link is that of the functions of the family, so an object fitted by an
# earlier version under a renamed link is drawn on the scale of its fits.
.trend_logs <- function(scmpObj, logs) {
    if (is.null(logs)) {
        logs <- .log_link_fits(scmpObj)
    }
    logs
}

# TRUE when the family of the fits has the functions of the log link
.log_link_fits <- function(scmpObj) {
    identical(.true_link(scmpObj@Parameters@distribution), "log")
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

    # Divide the points by the size factors when the fits used the offset,
    # which only log-link fits can
    use_offset <- isTRUE(scmpObj@Parameters@offset) &&
        length(offset_vector) > 0L && .log_link_fits(scmpObj)
    if (use_offset) {
        # Why: the curve is the model without the offset, the mean at a
        # size factor of 1, so the points are scaled to that depth too.
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

    if (use_offset) {
        line.df <- points.df
    }

    list(
        points = points.df, lines = line.df, curves = curve.df,
        path.names = path.names
    )
}
