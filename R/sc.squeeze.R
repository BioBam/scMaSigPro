#' @title Pseudo-bulking with optimal number of Pseudotime based bins
#'
#' @description
#' `sc.squeeze()` discretizes a continuous time series column into bins
#' of equal size using entropy-based binning method. It automatically calculates
#' the optimal number of bins using one of the supported methods.
#'
#' @importFrom parallel mclapply
#' @importFrom entropy discretize
#' @importFrom dplyr bind_rows group_by_at
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom dplyr group_by filter ungroup summarise summarize
#' @importFrom SummarizedExperiment assay colData<-
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values in 'Sparse' data. (Default is "Pseudotime").
#' @param path_col A character string representing the column name for branching
#' path assignment in 'Sparse' or 'Dense' data. (Default is `path_prefix`).
#' @param bin_ptime_col A character string representing the column name
#' for binned Pseudotime values in 'Dense' data.
#' (Default is "scmp_binned_pseudotime").
#' @param bin_mem_col A character string representing the name of the column in
#' which cells per bin are stored. (Default is "scmp_bin_members").
#' @param bin_col A character string representing the name of the column in which
#' bin labels are stored. (Default is "scmp_bin").
#' @param bin_size_col A character string representing the name of the column in
#' which bin sizes per bin are stored. (Default is "scmp_bin_size").
#' @param bin_method  A character string representing the algorithm used for
#' binning. Available options: "Freedman.Diaconis",
#' "Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal". (Default = "Sturges")
#' @param drop_fac A numeric value specifying the factor by which to adjust the
#' number of bins if the initial binning results in too many/few bins.
#' (Default = 1).
#' @param assay_name Name of the Assay in sparse data from which the counts are
#' used. (Default = "counts").
#' @param split_bins If bin sizes are greater than mean + sd, split the bin into
#' smaller bins by re-running the sc.squeeze() function. (Default = FALSE).
#' @param prune_bins If bin sizes are smaller than mean - sd, remove the bin.
#' (Default = FALSE).
#' @param drop_trails If the paths have different lengths of the binned pseudotime,
#' drop the bins from the path with more bins. (Default = FALSE).
#' @param fill_gaps If corresponding bin is missing for a time-point, pull the
#' successive bins and fill the gaps.
#' @param aggregate A character string specifying the method to aggregate counts
#' within each cluster. Available options are 'mean' or 'sum'. (Default = "sum").
#' @param verbose  Print detailed output in the console. (Default is TRUE)
#' @param additional_params Pass additional parameters as a named list.
#' See examples
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with updated `Dense`
#' slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso \code{\link{estBinSize}}, \code{\link[entropy]{discretize}},
#' \code{\link{create_range}}
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
    assertthat::assert_that(is(scmpObj, "ScMaSigPro"),
        msg = "Please provide object of class 'ScMaSigPro'."
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
    assertthat::assert_that(
        all(
            assay_name %in% names(scmpObj@Sparse@assays@data@listData)
        ),
        msg = paste0("'", assay_name, "' ", "doesn't exist in scmpObj.")
    )

    # Checks
    assertthat::assert_that(ptime_col %in% colnames(raw_cell_metadata),
        msg = paste0(
            "'", ptime_col, "' does not exist in cell.level.metadata. ",
            "Please review the 'ptime_col' parameter."
        )
    )
    assertthat::assert_that(path_col %in% colnames(raw_cell_metadata),
        msg = paste0(
            "'", path_col,
            "' does not exist in cell.level.metadata. Please review the 'path_col' parameter."
        )
    )
    assertthat::assert_that(drop_fac >= 0.3,
        msg = "Invalid value for 'drop_fac'. It should be between 0.3 and 1."
    )
    assertthat::assert_that(
        all(
            bin_method %in% c(
                "Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane", "Scott.Normal"
            )
        ),
        msg = "Available binning methods are 'Freedman.Diaconis', 'Sqrt', 'Sturges', 'Rice', 'Doane', and 'Scott.Normal'"
    )
    if (!is.null(additional_params)) {
        assertthat::assert_that(is.list(additional_params),
            msg = paste(
                "Please provide 'additional_params' as a named list.",
                "See details for more information"
            )
        )

        assertthat::assert_that(names(additional_params) %in% c("use_unique_time_points"),
            msg = "Allowed additional parameters are 'use_unique_time_points'."
        )
    }

    # Add a column
    raw_cell_metadata$cell <- rownames(raw_cell_metadata)

    # Get the avaible paths
    avail.paths <- as.vector(unique(raw_cell_metadata[[path_col]]))

    # Check for path
    assertthat::assert_that(length(avail.paths) >= 2,
        msg = "Invalid number of paths detected. Please make sure that dataset has at least two paths"
    )

    if (verbose) {
        message("Computing optimal bin-size with ", bin_method, " method.")
        message(
            "Number of available path in the dataset: ", length(avail.paths)
        )
        message("Paths: ", paste(avail.paths, collapse = ", "))
        message("Drop factor: ", drop_fac)
    }

    # Order
    avail.paths <- avail.paths[order(avail.paths)]
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
                                                  bin.members.colname = bin_mem_col,
                                                  v = verbose,
                                                  use.unique.time.points = additional_params$use_unique_time_points,
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
                message("Error message: ", e$message)
                stop("Unable to estimate bin size")
            }
        )

        # Get Bin Table
        bin_table <- extract_interval(
            time.vector = time_vector,
            nBins = estBins,
            bin = bin, bin.size = bin.size, lbound = lbound,
            ubound = ubound
        )

        # Client-Verbose
        if (verbose) {
            message(
                "For ", path, ", ", length_n,
                " time points has been compressed to ", nrow(bin_table), " bins"
            )
        }

        # Get Mean and SD
        mean_value <- round(mean(bin_table[[bin.size]]))
        sd_value <- round(sd(bin_table[[bin.size]]))

        # Get thresholds
        max.allowed <- abs(mean_value + sd_value)

        # Initate Runners
        new_max <- max.allowed + 1
        it <- 1

        # Split bins
        if (split) {
            if (verbose) {
                message(
                    "Optimizing bin sizes, with maximum allowed bin size as ",
                    max.allowed
                )
            }

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
                    bin = bin,
                    drop = drop_factor,
                    bin.size = bin.size,
                    method = method.bin
                )

                # Update new_max and new_min after the optimization step
                new_max <- as.numeric(max(bin_table[[bin.size]]))

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

        # Loop over each row of raw_cell_metadata
        for (i in seq_len(nrow(path.frame))) {
            # Get row from the Path frame
            raw_i <- path.frame[i, , drop = FALSE]
            pseudotime_i <- raw_i[[ptime_col]]

            # Check limits
            in_range <- pseudotime_i >= min(bin_table[, lbound]) &&
                pseudotime_i <= max(bin_table[, ubound])
            if (in_range) {
                # Check if the last elemenst (Last interval close both)
                if (max(path.frame[[ptime_col]] == pseudotime_i)) {
                    matching_bins <- bin_table[
                        (
                            pseudotime_i >= bin_table[, lbound] &
                                pseudotime_i <= bin_table[, ubound]
                        ), ,
                        drop = FALSE
                    ]
                } else {
                    # Filter bin_table based on pseudotime_value (right open left close)
                    matching_bins <- bin_table[
                        (
                            pseudotime_i >= bin_table[, lbound] &
                                pseudotime_i < bin_table[, ubound]
                        ), ,
                        drop = FALSE
                    ]
                }
            } else if (pseudotime_i > max(bin_table[, ubound])) {
                matching_bins <- bin_table[bin_table[, ubound] == max(bin_table[, ubound]), , drop = FALSE]
            } else if (pseudotime_i < min(bin_table[, lbound])) {
                matching_bins <- bin_table[bin_table[, lbound] == min(bin_table[, lbound]), , drop = FALSE]
            }

            if (nrow(matching_bins) > 0) {
                # If two match present take the first one
                if (nrow(matching_bins) == 2) {
                    path.frame[i, bin_ptime_col] <- as.numeric(matching_bins[[bin_ptime_col]][1])
                } else {
                    path.frame[i, bin_ptime_col] <- as.numeric(matching_bins[[bin_ptime_col]])
                }
            }
        }

        # Convert result to a data frame
        processed.path.frame <- as.data.frame(path.frame)
        processed.path.frame <- processed.path.frame[!is.na(processed.path.frame[[bin_ptime_col]]), , drop = FALSE]

        # Set the 'cell' column as rownames
        rownames(processed.path.frame) <- processed.path.frame$cell

        # Create Binned Path Frame
        binned.path.frame <- processed.path.frame %>%
            group_by_at(bin.time.col) %>%
            summarise(!!sym(bin.members.colname) := paste0(!!sym(cell), collapse = "|"))

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

    # List to frame
    processed_cell_metadata <- bind_rows(processed_cell_metadata.list) %>%
        as.data.frame()
    processed_binned_cell_metadata <- bind_rows(processed_binned_cell_metadata.list) %>%
        as.data.frame()

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

    # Prune and Trails
    if (prune_bins) {
        mean_bin_size <- mean(processed_binned_cell_metadata[[bin_size_col]])
        sd_bin_size <- sd(processed_binned_cell_metadata[[bin_size_col]])
        min_size_allowed <- abs(mean_bin_size - sd_bin_size)
        processed_binned_cell_metadata <- processed_binned_cell_metadata[processed_binned_cell_metadata[[bin_size_col]] >= min_size_allowed, , drop = FALSE]
    }

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
        pB.frame_tmp <- pB.frame %>%
            group_by(!!sym(path_col)) %>%
            summarize(max_time = max(!!sym(bin_ptime_col))) %>%
            ungroup()

        # Find the minimum of the max 'scmp_binned_pseudotime' across all 'Path'
        min_max_time <- min(pB.frame_tmp$max_time)

        # Identify rows to be removed
        rows_to_remove <- pB.frame %>%
            filter(!!sym(bin_ptime_col) > min_max_time)

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
        pB.frame <- pB.frame %>%
            filter(!!sym(bin_ptime_col) <= min_max_time)

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
    methods::validObject(scmpObj)
    return(scmpObj)
}
