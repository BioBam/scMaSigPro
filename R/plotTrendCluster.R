#' @title Plot the trends of the clusters
#'
#' @description
#' `plotTrendCluster()` plots, in one panel per cluster of
#' [sc.cluster.trend()], the points of its genes along the binned pseudotime
#' and a loess curve per path through their fitted curves.
#'
#' @param ylab The label of the y axis. The default, `NULL`, is
#'   `"log(Pseudobulk Expression)"` for log-link fits and when the points are
#'   logged, and `"Pseudobulk Expression"` otherwise.
#' @param logs Whether to log the points, as in [plotTrend()]. The default,
#'   `NULL`, logs them when the fits used the log link.
#' @param significant Whether the genes must be in a gene list of
#'   [sc.filter()]; the others are left out with a warning.
#' @param loess_span The span of the loess curves.
#'
#' @importFrom ggplot2 .data aes element_blank element_line element_text
#'   facet_wrap geom_path geom_point geom_smooth ggplot labs scale_color_manual
#'   theme theme_classic xlab ylab
#'
#' @return A ggplot object with one panel per cluster.
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @family plots
#' @inheritParams scmp-params
#' @inheritParams plotTrend
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' plotTrendCluster(scmp.ob, verbose = FALSE)
#'
#' @export
plotTrendCluster <- function(scmpObj,
                             xlab = "Pooled Pseudotime",
                             ylab = NULL,
                             summary_mode = "median",
                             logs = NULL, logType = "log",
                             smoothness = 1,
                             verbose = TRUE,
                             pseudoCount = 1,
                             significant = FALSE,
                             curves = TRUE,
                             lines = FALSE,
                             points = TRUE,
                             loess_span = 0.8) {
    # Global vars
    scmp_clusters <- "scmp_clusters"
    feature_id <- "feature_id"
    # Check Object Validity
    .check_scmp(scmpObj)
    logs <- .trend_logs(scmpObj, logs)
    if (is.null(ylab)) {
        ylab <- if (logs || .log_link_fits(scmpObj)) {
            "log(Pseudobulk Expression)"
        } else {
            "Pseudobulk Expression"
        }
    }

    # Check
    .check(!isEmpty(scmpObj@Significant@clusters),
        msg = "Please run 'sc.cluster.trend', before plotting cluster trends"
    )

    # Check Assertion
    .check(curves || lines || points,
        msg = "At least one of 'curves', 'lines', or 'points' must be TRUE."
    )

    summary_mode <- .check_choice(summary_mode, c("median", "mean"),
        msg = paste(
            paste0("'", summary_mode, "'"),
            "is not a valid option. Please use one of",
            paste(c("median", "mean"), collapse = ", ")
        )
    )
    if (logs) {
        logType <- .check_choice(logType, c("log", "log2", "log10"),
            msg = "'logType' should be one of 'log2', 'log10', 'log'"
        )
    }

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

    # The layers of plotTrend() for each gene, with the columns renamed by
    # name: x_axis, y_axis and group (the path)
    axes <- function(df, x, y) {
        names(df)[match(c(x, y, "path"), names(df))] <- c(
            "x_axis", "y_axis", "group"
        )
        df
    }
    genes <- gene.cluster.df[[feature_id]]
    data.list <- lapply(genes, function(gene_i) {
        trend <- tryCatch(
            .trendData(scmpObj, gene_i,
                smoothness = smoothness, logs = logs, logType = logType,
                pseudoCount = pseudoCount, significant = significant,
                summary_mode = summary_mode
            ),
            error = function(e) NULL
        )
        if (is.null(trend)) {
            return(NULL)
        }
        list(
            points = if (points) axes(trend$points, "pooled.time", "pb.counts"),
            line = if (lines) axes(trend$lines, "pooled.time", "pb.counts"),
            curve = if (curves) axes(trend$curves, "x", "y")
        )
    })
    names(data.list) <- genes

    if (verbose) {
        message("Computed trends for each gene, collapsing trends...")
    }

    # Genes without a trend, such as genes that are not significant with
    # significant = TRUE, are left out with a warning
    failed <- genes[vapply(data.list, is.null, logical(1))]
    if (length(failed) > 0) {
        warning(
            "No trend for ", length(failed), " gene(s), left out of the ",
            "plot: ", paste(failed, collapse = ", ")
        )
    }
    data.list <- Filter(Negate(is.null), data.list)

    # Create List as per the clusters, keyed by cluster
    cluster_ids <- unique(gene.cluster.df[[scmp_clusters]])
    cluster.data.list <- lapply(cluster_ids, function(clus) {
        # Get gene names with same cluster
        in_cluster <- gene.cluster.df[[scmp_clusters]] == clus
        gene_set <- gene.cluster.df[in_cluster, , drop = FALSE][[feature_id]]

        # Subset the datalist
        data.list[gene_set]
    })
    names(cluster.data.list) <- cluster_ids

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


        # Summarise y_axis within each x_axis and group (path)
        summarise_y <- function(df, fun) {
            .summarise_by(df, c("x_axis", "group"), "y_axis", fun, na.rm = TRUE)
        }
        if (summary == "mean") {
            if (!is.null(line.df)) {
                line.df <- summarise_y(line.df, mean)
            }
            if (!is.null(point.df)) {
                point.df <- summarise_y(point.df, mean)
            }
            if (!is.null(curve.df)) {
                curve.df <- summarise_y(curve.df, mean)
            }
        } else if (summary == "median") {
            if (!is.null(line.df)) {
                line.df <- summarise_y(line.df, median)
            }
            if (!is.null(point.df)) {
                point.df <- summarise_y(point.df, median)
            }
            if (!is.null(curve.df)) {
                curve.df <- summarise_y(curve.df, median)
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
        cluster_name <- names(collapsed.df)[i]

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

    # Cluster number: the second "_"-separated piece of "cluster_<n>"; the
    # column is a factor after merge(), hence as.character()
    cluster_no <- function(x) {
        vapply(strsplit(as.character(x), "_", fixed = TRUE), "[", "", 2)
    }

    # Combine the data frames of the same type
    if (lines) {
        lines_combined <- as.data.frame(do.call(rbind, lines_list))
        lines_combined <- merge(freq.table, lines_combined, by = "cluster")
        lines_combined[["cluster"]] <- paste0(
            paste("Cluster", cluster_no(lines_combined[["cluster"]]),
                sep = ": "
            ),
            " (", lines_combined[["num"]], " Features)"
        )
    }
    if (curves) {
        curves_combined <- as.data.frame(do.call(rbind, curves_list))
        curves_combined <- merge(freq.table, curves_combined, by = "cluster")
        curves_combined[["cluster"]] <- paste0(
            paste("Cluster", cluster_no(curves_combined[["cluster"]]),
                sep = ": "
            ),
            " (", curves_combined[["num"]], " Features)"
        )
    }
    if (points) {
        points_combined <- as.data.frame(do.call(rbind, points_list))
        points_combined <- merge(freq.table, points_combined, by = "cluster")
        # Fix cluster info
        points_combined[["cluster"]] <- paste0(
            paste("Cluster", cluster_no(points_combined[["cluster"]]),
                sep = ": "
            ),
            " (", points_combined[["num"]], " Features)"
        )
    }

    if (verbose) {
        message("Calculating 'loess', hang tight this can take a while depending on the smoothness...")
    }

    # Initiate plotting
    p <- ggplot()

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
            formula = y ~ x,
            method = "loess", span = loess_span,
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
    return(p)
}
