#' @title Bin the cells of each path along pseudotime
#'
#' @description
#' `sc.squeeze()` divides the pseudotime range of each path into bins of
#' equal width and adds up the counts of the cells of each bin into
#' pseudo-bulk counts. The number of bins of a path is that of a histogram
#' rule ([estBinSize()]) times `drop_fac`, rounded up.
#'
#' Each cell goes to the bin whose interval holds its pseudotime. As in
#' Srivastava et al. (2024), Section 2.1, "bins are structured as left-closed
#' and right-open intervals, except for the last bin, which is closed on both
#' ends": a cell on the bound of two bins goes to the later one. The bounds
#' are stored in the `scmp_l_bound` and `scmp_u_bound` columns of the bin
#' metadata ([cDense()]).
#'
#' The paths are binned, and their bins stored, in the order of the levels of
#' the path column when it is a factor, and otherwise in the C-locale (byte)
#' order of the labels, whatever the collation of the session.
#' [sc.set.poly()] takes the first path as the reference unless told
#' otherwise.
#'
#' @param bin_method The rule for the number of bins: `"Sturges"`, `"Sqrt"`,
#'   `"Rice"`, `"Doane"`, `"Freedman.Diaconis"` or `"Scott.Normal"` (see
#'   [estBinSize()]).
#' @param drop_fac The drop factor of Srivastava et al. (2024), Section 2.1,
#'   at least 0.3: it multiplies the number of bins of the rule, and the
#'   product is rounded up. Values below 1 give fewer, wider bins.
#' @param split_bins Whether to divide the bins of a path that hold more cells
#'   than the mean plus the standard deviation of its bin sizes (see Details).
#' @param prune_bins Whether to remove the bins of a path that hold fewer
#'   cells than the mean minus the standard deviation of its bin sizes (see
#'   Details).
#' @param drop_trails Whether to drop the last bins of the longer paths, so
#'   that every path has the binned pseudotime of the shortest.
#' @param fill_gaps Whether to renumber the bins of a path whose binned
#'   pseudotime has gaps, such as after pruning.
#' @param additional_params A list: `use_unique_time_points = TRUE` bins the
#'   distinct pseudotime values of each path instead of its cells.
#'
#' @details
#' `split_bins` and `prune_bins` follow Section 2.2 of Srivastava et al.
#' (2024): "scMaSigPro calculates the mean (\eqn{\mu}) and standard deviation
#' (\eqn{\sigma}) of the sizes of all bins per path to establish maximum
#' (\eqn{\mu+\sigma}) and minimum (\eqn{\mu-\sigma}) allowable bin sizes. Bins
#' larger than the maximum size are divided, and those smaller than the
#' minimum size can be optionally removed, ensuring uniform bin sizes. The
#' process is repeated until all bins conform to the established size
#' criteria." Both limits are set once, from the bins of the rule, empty bins
#' included. A bin above the maximum is divided into the number of bins that
#' the rule and `drop_fac` give for its cells, until no bin is larger;
#' `sc.squeeze()` stops when a division gives back bins it has already made.
#' A path with a single bin keeps it. Pruning comes after the divisions: a
#' negative minimum removes no bin, and the cells of a removed bin stay in
#' the `Sparse` slot without a bin.
#'
#' @importFrom SummarizedExperiment assay colData<-
#'
#' @return The object with the bins and their pseudo-bulk counts
#' (`bulk.counts`) in the `Dense` slot and the binned pseudotime of each cell
#' in the cell metadata.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso [estBinSize()] for the rules that set the number of bins;
#' [plotBinTile()] to see the bins.
#'
#' @family workflow steps
#' @inheritParams scmp-params
#'
#' @references
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' @examples
#' data("splat.sim", package = "scMaSigPro")
#' scmp <- as_scmp(splat.sim,
#'     additional_params = list(
#'         labels_exist = TRUE, exist_ptime_col = "Step",
#'         exist_path_col = "Group"
#'     ),
#'     verbose = FALSE
#' )
#' scmp <- sc.squeeze(scmp)
#' bins <- cDense(scmp)
#' head(bins[, c("Path", "scmp_binned_pseudotime", "scmp_bin_size")])
#'
#' @export
sc.squeeze <- function(scmpObj,
                       ptime_col = scmpObj@Parameters@ptime_col,
                       path_col = scmpObj@Parameters@path_col,
                       bin_method = "Sturges",
                       drop_fac = 1,
                       verbose = FALSE,
                       bin_mem_col = "scmp_bin_members",
                       bin_col = "scmp_bin",
                       bin_size_col = "scmp_bin_size",
                       bin_ptime_col = "scmp_binned_pseudotime",
                       split_bins = FALSE,
                       prune_bins = FALSE,
                       assay_name = "counts",
                       drop_trails = FALSE,
                       aggregate = "sum",
                       fill_gaps = FALSE,
                       additional_params = list(use_unique_time_points = FALSE)) {
    # Initiate Variable
    scmp_bin_lower_bound <- "scmp_l_bound"
    scmp_bin_upper_bound <- "scmp_u_bound"
    cell <- "cell"

    # Check Object Validity
    .check_scmp(scmpObj)
    .check_strings(
        ptime_col = ptime_col, path_col = path_col, bin_mem_col = bin_mem_col,
        bin_col = bin_col, bin_size_col = bin_size_col,
        bin_ptime_col = bin_ptime_col, assay_name = assay_name
    )

    # Extract cell metadata
    raw_cell_metadata <- as.data.frame(colData(scmpObj@Sparse))

    # Drop Columns if exist
    cols_to_drop <- c(
        scmpObj@Parameters@bin_size_col,
        scmpObj@Parameters@bin_ptime_col,
        "scmp_u_bound", "scmp_l_bound"
    )
    raw_cell_metadata <- raw_cell_metadata[, !colnames(
        raw_cell_metadata
    ) %in% cols_to_drop, drop = FALSE]

    # Count slot
    .check(
        all(
            assay_name %in% assayNames(scmpObj@Sparse)
        ),
        msg = paste0("'", assay_name, "' ", "doesn't exist in scmpObj.")
    )

    # Checks
    .check(ptime_col %in% colnames(raw_cell_metadata),
        msg = paste0(
            "'", ptime_col, "' does not exist in cell.level.metadata. ",
            "Please review the 'ptime_col' parameter."
        )
    )
    .check(path_col %in% colnames(raw_cell_metadata),
        msg = paste0(
            "'", path_col, "' does not exist in cell.level.metadata. ",
            "Please review the 'path_col' parameter."
        )
    )
    # Why: Section 2.1 of the paper bounds the drop factor below by 0.3.
    .check(drop_fac >= 0.3,
        msg = "Invalid value for 'drop_fac'. It should be at least 0.3."
    )
    bin_method <- .check_choice(bin_method,
        c(
            "Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane",
            "Scott.Normal"
        ),
        msg = paste(
            "Available binning methods are 'Freedman.Diaconis', 'Sqrt',",
            "'Sturges', 'Rice', 'Doane', and 'Scott.Normal'"
        )
    )
    aggregate <- .check_choice(aggregate, c("sum", "mean"),
        msg = "Invalid 'aggregate' value. Please choose either 'mean' or 'sum'."
    )
    # NULL or an empty list keep the default
    if (is.null(additional_params)) {
        additional_params <- list()
    }
    .check(is.list(additional_params),
        msg = paste(
            "Please provide 'additional_params' as a named list.",
            "See details for more information"
        )
    )
    .check(
        length(additional_params) == 0L ||
            identical(names(additional_params), "use_unique_time_points"),
        msg = "Allowed additional parameters are 'use_unique_time_points'."
    )
    use_unique <- additional_params[["use_unique_time_points"]]
    if (is.null(use_unique)) {
        use_unique <- FALSE
    }

    # Add a column
    raw_cell_metadata$cell <- rownames(raw_cell_metadata)

    # Get the avaible paths
    avail.paths <- as.vector(unique(raw_cell_metadata[[path_col]]))

    # Check for path
    .check(length(avail.paths) >= 2,
        msg = paste(
            "Invalid number of paths detected. Please make sure that dataset",
            "has at least two paths"
        )
    )

    if (verbose) {
        message("Computing optimal bin-size with ", bin_method, " method.")
        message(
            "Number of available path in the dataset: ", length(avail.paths)
        )
        message("Paths: ", paste(avail.paths, collapse = ", "))
        message("Drop factor: ", drop_fac)
    }

    # Order the paths: by the levels of a factor, as model formulas do, and
    # otherwise by the bytes of the labels (the C locale) in every locale:
    # order() alone follows the session's collation, which put "root" first
    # in en_US and after "Y_15" in C
    path_values <- raw_cell_metadata[[path_col]]
    if (is.factor(path_values)) {
        avail.paths <- intersect(levels(path_values), avail.paths)
    } else {
        avail.paths <- avail.paths[order(avail.paths, method = "radix")]
    }
    if (verbose) {
        message("Setting ", avail.paths[1], " as reference.")
    }

    # Apply transformations on data
    discrete.list <- lapply(avail.paths, function(path,
                                                  design.frame = raw_cell_metadata,
                                                  drop_factor = drop_fac,
                                                  path.col = path_col,
                                                  bin.size = bin_size_col,
                                                  bin = bin_col,
                                                  time.col = ptime_col,
                                                  method.bin = bin_method,
                                                  bin.time.col = bin_ptime_col,
                                                  split = split_bins,
                                                  prune = prune_bins,
                                                  bin.members.colname = bin_mem_col,
                                                  v = verbose,
                                                  use.unique.time.points = use_unique,
                                                  lbound = scmp_bin_lower_bound,
                                                  ubound = scmp_bin_upper_bound) {
        # Get the cells belonging to path
        path.frame <- design.frame[design.frame[[path.col]] == path, , drop = FALSE]

        # Extract the time information as a vector
        time_vector <- path.frame[, time.col]
        length_n <- length(time_vector)


        if (use.unique.time.points) {
            time_vector <- unique(time_vector)
            length_n <- length(time_vector)
            if (v) {
                message("Using only unique points in the time series")
            }
        }

        # Validation
        if (length_n <= 7) {
            message("Time points are already less than 7 in ", path)
        }

        # Calculate Optimal Number of Bins
        tryCatch(
            expr = {
                # Why: Supplementary Section S1.1 of the paper rounds the
                # rule times the drop factor up.
                estBins <- ceiling(estBinSize(
                    time_vector = time_vector, nPoints = length_n,
                    drop_fac = drop_factor, bin_method = method.bin
                ))

                if (verbose) {
                    message(
                        "Estimated Bin Sizes = ", estBins, " with ", bin_method,
                        " binning for ", length_n, " time points for ", path
                    )
                }
            },
            error = function(e) {
                message("Could not estimate the number of bins: ", e$message)
                stop("Unable to estimate bin size")
            }
        )

        # Get Bin Table
        bin_table <- extract_interval(
            time.vector = time_vector,
            nBins = estBins,
            bin.size = bin.size, lbound = lbound, ubound = ubound
        )

        # Client-Verbose
        if (verbose) {
            message(
                "For ", path, ", ", length_n,
                " time points has been compressed to ", nrow(bin_table), " bins"
            )
        }

        # Section 2.2: the mean and the standard deviation of the sizes of
        # all bins of the path, empty intervals included, set the largest
        # and the smallest sizes allowed, mean + sd and mean - sd, once,
        # before any split
        mean_value <- mean(bin_table[[bin.size]])
        sd_value <- sd(bin_table[[bin.size]])
        max.allowed <- mean_value + sd_value
        min.allowed <- mean_value - sd_value

        # Initate Runners
        new_max <- max.allowed + 1
        it <- 1

        # Split bins; the sizes of a path with a single bin have no standard
        # deviation, and its bin is kept
        if (split && nrow(bin_table) > 1) {
            if (verbose) {
                message(
                    "Optimizing bin sizes, with maximum allowed bin size as ",
                    round(max.allowed, 2)
                )
            }

            # optimize_bin_max() is deterministic: a bin table seen before
            # means that the loop would never end
            seen <- list(bin_table)

            # Adjust maximum Size
            while (new_max > max.allowed) {
                if (verbose) {
                    message("Iteration ", it)
                }

                # Call your 'optimize.bin.width' function
                bin_table <- optimize_bin_max(
                    bin_table = bin_table,
                    max_allowed = max.allowed,
                    verbose = verbose,
                    time_vector = time_vector,
                    lbound = lbound,
                    ubound = ubound,
                    drop = drop_factor,
                    bin.size = bin.size,
                    method = method.bin
                )

                # Update new_max and new_min after the optimization step
                new_max <- as.numeric(max(bin_table[[bin.size]]))
                if (new_max > max.allowed) {
                    repeated <- vapply(seen, identical, logical(1), bin_table)
                    .check(!any(repeated),
                        msg = paste0(
                            "Splitting the bins of '", path, "' with ",
                            "bin_method = '", method.bin, "' makes no ",
                            "progress: the largest bin keeps ", new_max,
                            " cells, more than the ", round(max.allowed, 2),
                            " allowed. Use split_bins = FALSE or another ",
                            "bin_method."
                        )
                    )
                    seen <- c(seen, list(bin_table))
                }

                # Increment the iteration counter
                it <- it + 1
            }
            if (verbose) {
                message(
                    "Finally, for ", path, ", ", length_n,
                    " time points has been compressed to ", nrow(bin_table),
                    " bins and the sum is ", sum(bin_table[[bin.size]])
                )
            }
        }

        rownames(bin_table) <- NULL
        bin_table[[bin_ptime_col]] <- as.numeric(rownames(bin_table))

        # The bin of every cell, by the numeric bounds; a cell that no bin
        # holds keeps NA and is removed below
        cell_bin <- .assign_bins(
            path.frame[[ptime_col]], bin_table[[lbound]], bin_table[[ubound]]
        )
        path.frame[[bin_ptime_col]] <- bin_table[[bin_ptime_col]][cell_bin]

        # Convert result to a data frame
        processed.path.frame <- as.data.frame(path.frame)
        processed.path.frame <- processed.path.frame[!is.na(processed.path.frame[[bin_ptime_col]]), , drop = FALSE]

        # Set the 'cell' column as rownames
        rownames(processed.path.frame) <- processed.path.frame$cell

        # Create Binned Path Frame: one row per binned time, in increasing
        # order, with the cells of the bin (in metadata order) joined by "|"
        cell_times <- processed.path.frame[[bin.time.col]]
        bin_times <- sort(unique(cell_times))
        bin_members <- vapply(
            split(processed.path.frame[[cell]], match(cell_times, bin_times)),
            paste0, character(1),
            collapse = "|", USE.NAMES = FALSE
        )
        binned.path.frame <- data.frame(bin_times, bin_members)
        colnames(binned.path.frame) <- c(bin.time.col, bin.members.colname)

        # Missing columns
        tmp.bin.name <- paste0(path, "_bin_", seq(1, nrow(binned.path.frame)))
        binned.path.frame[[bin]] <- tmp.bin.name
        binned.path.frame[[path.col]] <- path


        tmp.bin.size <- apply(binned.path.frame, 1, calc_bin_size,
            clus_mem_col = bin.members.colname
        )
        binned.path.frame[[bin.size]] <- tmp.bin.size

        # Subset the bin tables
        bin_table_sub <- bin_table[, c(lbound, ubound, bin.time.col)]

        # Merge
        binned.path.frame <- merge(binned.path.frame, bin_table_sub, by = bin.time.col)

        # Prune the bins smaller than mean - sd: a negative minimum removes
        # none, and a path with a single bin, whose size has no sd, keeps it
        if (prune && !is.na(min.allowed)) {
            binned.path.frame <- binned.path.frame[
                binned.path.frame[[bin.size]] >= min.allowed, ,
                drop = FALSE
            ]
        }

        return(list(
            sparse = processed.path.frame,
            compressed = binned.path.frame
        ))
    })

    # Bind rows and convert to data frame, then drop 'cell' column
    processed_cell_metadata.list <- lapply(discrete.list, function(element) {
        return(element[["sparse"]])
    })
    processed_binned_cell_metadata.list <- lapply(discrete.list, function(element) {
        return(element[["compressed"]])
    })

    # List to frame: the per-path frames have the same columns; the row
    # names are set from the 'cell' and bin columns below
    processed_cell_metadata <- do.call(rbind, processed_cell_metadata.list)
    processed_binned_cell_metadata <- do.call(
        rbind, processed_binned_cell_metadata.list
    )

    # Put the cell metadata back in the column order of the counts: the
    # per-path frames above were stacked path by path
    rownames(processed_cell_metadata) <- processed_cell_metadata$cell
    assay_cells <- colnames(assay(scmpObj@Sparse, assay_name, withDimnames = FALSE))
    if (is.null(assay_cells)) {
        assay_cells <- colnames(scmpObj@Sparse)
    }
    is_binned <- assay_cells %in% rownames(processed_cell_metadata)
    if (!all(is_binned)) {
        warning(
            sum(!is_binned), " cell(s) could not be assigned to a bin and are ",
            "removed from the 'Sparse' slot."
        )
        scmpObj@Sparse <- scmpObj@Sparse[, is_binned]
        assay_cells <- assay_cells[is_binned]
    }
    colData(scmpObj@Sparse) <- DataFrame(processed_cell_metadata[
        match(assay_cells, rownames(processed_cell_metadata)), ,
        drop = FALSE
    ])

    rownames(processed_binned_cell_metadata) <- processed_binned_cell_metadata[[bin_col]]

    if (fill_gaps) {
        processed_binned_cell_metadata_tmp <- data.frame()

        for (path in unique(processed_binned_cell_metadata[[path_col]])) {
            # Get path
            binned.path.frame <- processed_binned_cell_metadata[processed_binned_cell_metadata[[path_col]] == path, , drop = FALSE]

            # Generate refernce
            ref_time <- seq_len(nrow(binned.path.frame))

            # Align with the binned Pseudotime
            binned.pseudotime <- binned.path.frame[[bin_ptime_col]]

            # Check
            if (!all(ref_time == binned.pseudotime)) {
                binned.path.frame[[bin_ptime_col]] <- ref_time
                tmp.bin.name <- paste0(path, "_bin_", ref_time)
                binned.path.frame[[bin_col]] <- tmp.bin.name
                rownames(binned.path.frame) <- tmp.bin.name
            }

            processed_binned_cell_metadata_tmp <- rbind(
                processed_binned_cell_metadata_tmp,
                binned.path.frame
            )
        }
        processed_binned_cell_metadata <- processed_binned_cell_metadata_tmp
    }

    if (drop_trails) {
        pB.frame <- processed_binned_cell_metadata
        # Find the maximum 'scmp_binned_pseudotime' for each 'Path'
        max_time <- tapply(pB.frame[[bin_ptime_col]], pB.frame[[path_col]], max)

        # Find the minimum of the max 'scmp_binned_pseudotime' across all 'Path'
        min_max_time <- min(max_time)

        # Identify rows to be removed
        rows_to_remove <- pB.frame[
            which(pB.frame[[bin_ptime_col]] > min_max_time), ,
            drop = FALSE
        ]

        # Print a message about the rows that will be removed
        if (verbose) {
            if (nrow(rows_to_remove) > 0) {
                message(sprintf(
                    "Dropped trailing bin %s from %s\n",
                    rows_to_remove[[bin_ptime_col]], rows_to_remove[[path_col]]
                ), appendLF = FALSE)
            }
        }

        # Filter out rows where 'scmp_binned_pseudotime' is greater than 'min_max_time'
        pB.frame <- pB.frame[
            which(pB.frame[[bin_ptime_col]] <= min_max_time), ,
            drop = FALSE
        ]

        processed_binned_cell_metadata <- pB.frame
    }

    # Bin metadata only; pb_counts() adds the assay
    scmpObj@Dense <- SingleCellExperiment(
        colData = DataFrame(processed_binned_cell_metadata)
    )

    # Get Counts
    scmpObj <- pb_counts(
        scmpObj = scmpObj,
        bin_mem_col = bin_mem_col,
        bin_col = bin_col,
        assay_name = assay_name,
        cluster_count_by = aggregate
    )

    # Update Slots
    scmpObj@Parameters@ptime_col <- ptime_col
    scmpObj@Parameters@path_col <- path_col
    scmpObj@Parameters@bin_method <- bin_method
    scmpObj@Parameters@bin_ptime_col <- bin_ptime_col
    scmpObj@Parameters@bin_col <- bin_col
    scmpObj@Parameters@bin_mem_col <- bin_mem_col
    scmpObj@Parameters@bin_size_col <- bin_size_col
    validObject(scmpObj)
    return(scmpObj)
}
