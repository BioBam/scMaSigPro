#' @title Pseudo-bulking with optimal number of pseudotime based bins
#'
#' @description
#' `squeeze()` discretizes a continuous time series column into bins
#' of equal size using entropy-based binning method. It automatically calculates
#' the optimal number of bins using one of the supported methods. The bin sizes
#' are also calculated and merged with the input cell_metadata.
#'
#' @param scmpObject object of Class scMaSigPro. See \code{\link{scmp}}
#' for more details.
#' @param pseudotime_colname Name of the column in `cell.metadata` storing
#' Pseudotime values. Generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is "Pseudotime").
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path. Generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is `path_prefix`).
#' @param bin_pseudotime_colname Name of the column to store the computed Pseudotime
#' bins.
#' @param bin_members_colname Name of the column in the 'annotated_cell_metadata'
#' @param bin_colname Name of the bin column name
#' @param bin_size_colname Setting the name of the bin size column.
#' @param bin_method A character string specifying the method to use in order to
#' estimate the optimal number of bins. Available options: "Freedman.Diaconis",
#' "Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal". See \code{\link{estBinSize}}
#' for more details. (Default = "Sturges").
#' @param drop_fac A numeric value specifying the factor by which to decrease the
#' number of bins if the initial binning results in too many bins. (Default = 1).
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' per path iteratively. Options: "universal", "individual. (Default = "universal").
#' @param additional_params Pass additional parameters as a named list. See Details.
#' @param assay_name Name of the Assay in the assay_name object from which retrieve the counts.
#' (Default = "counts").
#' @param split_bins If bin sizes are greater than
#' @param prune_bins description
#' @param drop_trails description
#' @param fill_gaps description
#' @param cluster_count_by A character string specifying the method to use to
#' aggregate counts within each cluster. Available options are 'mean' or 'sum'. (Default = "sum").
#'
#' @return
#' A data.frame that contains the original data plus additional columns:
#' - 'bin' : Name of the bin
#' - 'bin_size' : Size of the bin
#' - 'binned_time' : Interval range of each bin
#' This function returns the merged data.frame with new discretized
#' pseudotime_colname, preserving the original rownames.
#'
#' @details
#' This function performs the following steps:
#' - Adds a new column 'cell' to the input data.frame which copies the row names.
#' - Extracts the time series data from the specified column of the input data.frame.
#' - Calculates the optimal number of bins using the specified method.
#' - Prints the estimated number of bins if verbose is set to TRUE.
#' - Discretizes the time series data into bins using the entropy-based binning method.
#' - Merges the original data.frame with the new binned time series data.
#' - Removes the 'cell' column and sets the row names back to the original row names of the input data.frame.
#' - Returns the merged data.frame.
#' @importFrom assertthat assert_that
#' @importFrom parallel mclapply detectCores
#' @importFrom entropy discretize
#' @importFrom dplyr left_join join_by mutate select bind_rows group_by_at summarise rename_with
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom dplyr group_by filter row_number ungroup summarise summarize
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso \code{\link{estBinSize}}, \code{\link{discretize}}, \code{\link{create_range}}
#'
#' @export

sc.squeeze <- function(scmpObject,
                       pseudotime_colname = scmpObject@param@pseudotime_colname,
                       path_colname = scmpObject@param@path_colname,
                       bin_method = "Sturges",
                       drop_fac = 1,
                       verbose = FALSE,
                       bin_members_colname = "scmp_bin_members",
                       bin_colname = "scmp_bin",
                       bin_size_colname = "scmp_bin_size",
                       bin_pseudotime_colname = "scmp_binned_pseudotime",
                       split_bins = FALSE,
                       prune_bins = FALSE,
                       assay_name = "counts",
                       drop_trails = FALSE,
                       cluster_count_by = "sum",
                       fill_gaps = FALSE,
                       additional_params = list(use_unique_time_points = FALSE)) {
  # Initiate Variable
  scmp_bin_lower_bound <- "scmp_l_bound"
  scmp_bin_upper_bound <- "scmp_u_bound"
  cell <- "cell"

  # Check Object Validity
  assert_that(is(scmpObject, "scmp"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Extract cell metadata
  raw_cell_metadata <- as.data.frame(colData(scmpObject@sparse))

  # Drop Columns if exist
  cols_to_drop <- c(
    scmpObject@param@bin_size_colname,
    scmpObject@param@bin_pseudotime_colname,
    "scmp_u_bound", "scmp_l_bound"
  )
  raw_cell_metadata <- raw_cell_metadata[, !colnames(raw_cell_metadata) %in% cols_to_drop, drop = FALSE]

  # Count slot
  assert_that(
    all(
      assay_name %in% names(scmpObject@sparse@assays@data@listData)
    ),
    msg = paste0("'", assay_name, "' ", "doesn't exit in scmpObject.")
  )

  # Checks
  assert_that(pseudotime_colname %in% colnames(raw_cell_metadata),
    msg = paste0("'", pseudotime_colname, "' does not exist in cell.level.metadata Please review the 'pseudotime_colname' parameter.")
  )
  assert_that(path_colname %in% colnames(raw_cell_metadata),
    msg = paste0("'", path_colname, "' does not exist in cell.level.metadata. Please review the 'path_colname' parameter.")
  )
  assert_that(drop_fac >= 0.3,
    msg = "Invalid value for 'drop_fac'. It should be between 0.3 and 1."
  )
  assert_that(
    all(
      bin_method %in% c("Freedman.Diaconis", "Sqrt", "Sturges", "Rice", "Doane", "Scott.Normal")
    ),
    msg = "Available binning methods are 'Freedman.Diaconis', 'Sqrt', 'Sturges', 'Rice', 'Doane', and 'Scott.Normal'"
  )
  if (!is.null(additional_params)) {
    assert_that(is.list(additional_params),
      msg = "Please provide 'additional_params' as a named list.
      See details for more information"
    )

    assert_that(names(additional_params) %in% c("use_unique_time_points"),
      msg = "Allowed additional parameters are 'use_unique_time_points'."
    )
  }

  # Add a column
  raw_cell_metadata$cell <- rownames(raw_cell_metadata)

  # Get the avaible paths
  avail.paths <- as.vector(unique(raw_cell_metadata[[path_colname]]))

  # Check for path
  assert_that(length(avail.paths) >= 2,
    msg = "Invalid number of paths detected. Please make sure that dataset has at least two paths"
  )

  if (verbose) {
    message(paste("Computing optimal bin-size with", bin_method, "method."))
    message(paste("Number of available path in the dataset:", length(avail.paths)))
    message(paste("Paths:", paste(avail.paths, collapse = ", ")))
    message(paste("Drop factor:", drop_fac))
  }

  # Order
  avail.paths <- avail.paths[order(avail.paths)]
  if (verbose) {
    message(paste("Setting", avail.paths[1], "as reference."))
  }

  # Apply transformations on data
  discrete.list <- lapply(avail.paths, function(path, design.frame = raw_cell_metadata,
                                                drop_factor = drop_fac, path.col = path_colname,
                                                bin.size = bin_size_colname, bin = bin_colname,
                                                time.col = pseudotime_colname, method.bin = bin_method,
                                                bin.time.col = bin_pseudotime_colname,
                                                split = split_bins, bin.members.colname = bin_members_colname,
                                                v = verbose, use.unique.time.points = additional_params$use_unique_time_points,
                                                lbound = scmp_bin_lower_bound, ubound = scmp_bin_upper_bound) {
    # Get the cells belonging to path
    path.frame <- design.frame[design.frame[[path.col]] == path, , drop = FALSE]

    # Extract the time information as a vector
    time_vector <- path.frame[, time.col]
    length_n <- length(time_vector)


    if (use.unique.time.points) {
      time_vector <- unique(time_vector)
      length_n <- length(time_vector)
      if (v) {
        message(paste("Using only unique points in the time series"))
      }
    }

    # Validation
    if (length_n <= 7) {
      message(paste("Time points are already less than 7 in", path))
    }

    # Calculate Optimal Number of Bins
    tryCatch(
      expr = {
        estBins <- ceiling(estBinSize(
          time_vector = time_vector, nPoints = length_n,
          drop_fac = drop_factor, bin_method = method.bin
        ))

        if (verbose) {
          message(paste(
            "Estimated Bin Sizes =", estBins, "with",
            bin_method, "binning for", length_n, "time points for", path
          ))
        }
      },
      error = function(e) {
        message(paste("Error message: ", e$message))
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
      message(paste(
        "For", path, ",", length_n, "time points has been compressed to", nrow(bin_table), "bins"
      ))
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
        message(paste("Optimizing bin sizes, with maximum allowed bin size as", max.allowed))
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
    }

    if (verbose) {
      message(paste("Optimizing bin sizes, with maximum allowed bin size as", max.allowed))
    }

    if (verbose) {
      message(paste(
        "Finally, for", path, ",", length_n, "time points has been compressed to", nrow(bin_table), "bins and the sum is ", sum(bin_table[[bin.size]])
      ))
    }

    rownames(bin_table) <- NULL
    bin_table[[bin_pseudotime_colname]] <- as.numeric(rownames(bin_table))

    # Loop over each row of raw_cell_metadata
    for (i in 1:nrow(path.frame)) {
      # Get row from the Path frame
      raw_i <- path.frame[i, , drop = FALSE]
      # print(paste("Completed raw_i <- path.frame[i, , drop = FALSE] for", path))


      # View(raw_i)
      # stop()
      pseudotime_i <- raw_i[[pseudotime_colname]]
      # print(paste("Completed pseudotime_i <- raw_i[[pseudotime_colname]] for", path))

      # Check limits
      if (pseudotime_i >= min(bin_table[, lbound]) & pseudotime_i <= max(bin_table[, ubound])) {
        # Check if the last elemenst (Last interval close both)
        if (max(path.frame[[pseudotime_colname]] == pseudotime_i)) {
          matching_bins <- bin_table[(
            pseudotime_i >= bin_table[, lbound] &
              pseudotime_i <= bin_table[, ubound]), , drop = FALSE]
        } else {
          # Filter bin_table based on pseudotime_value (right open left close)
          matching_bins <- bin_table[
            (
              pseudotime_i >= bin_table[, lbound] &
                pseudotime_i < bin_table[, ubound]), ,
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
          path.frame[i, bin_pseudotime_colname] <- as.numeric(matching_bins[[bin_pseudotime_colname]][1])
        } else {
          path.frame[i, bin_pseudotime_colname] <- as.numeric(matching_bins[[bin_pseudotime_colname]])
        }
      }
    }

    # Convert result to a data frame
    processed.path.frame <- as.data.frame(path.frame)
    processed.path.frame <- processed.path.frame[!is.na(processed.path.frame[[bin_pseudotime_colname]]), , drop = FALSE]

    # print(paste("processed.path.frame <- as.data.frame(path.frame)", path))

    # Set the 'cell' column as rownames
    rownames(processed.path.frame) <- processed.path.frame$cell

    # print(paste("rownames(processed.path.frame) <- processed.path.frame$cell", path))

    # print(paste("setted rownames for", path))

    # Create Binned Path Frame
    binned.path.frame <- processed.path.frame %>%
      group_by_at(bin.time.col) %>%
      summarise(!!sym(bin.members.colname) := paste0(!!sym(cell), collapse = "|"))

    # Missing columns
    tmp.bin.name <- paste0(path, "_bin_", seq(1, nrow(binned.path.frame)))
    binned.path.frame[[bin]] <- tmp.bin.name
    binned.path.frame[[path.col]] <- path


    tmp.bin.size <- apply(binned.path.frame, 1, calc_bin_size, clus_mem_col = bin.members.colname)
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

  ## Add Processed Cell Matadata back with slot update
  scmpObject@sparse@colData <- DataFrame(processed_cell_metadata)

  # Set the 'cell' column as rownames
  rownames(processed_cell_metadata) <- processed_cell_metadata$cell
  processed_cell_metadata <- processed_cell_metadata[,colnames(processed_cell_metadata) != "cell", drop=FALSE] 
  rownames(processed_binned_cell_metadata) <- processed_binned_cell_metadata[[bin_colname]]

  # Prune and Trails
  if (prune_bins) {
    mean_bin_size <- mean(processed_binned_cell_metadata[[bin_size_colname]])
    sd_bin_size <- sd(processed_binned_cell_metadata[[bin_size_colname]])
    min_size_allowed <- abs(mean_bin_size - sd_bin_size)
    processed_binned_cell_metadata <- processed_binned_cell_metadata[processed_binned_cell_metadata[[bin_size_colname]] >= min_size_allowed, , drop = FALSE]
  }

  if (fill_gaps) {
    processed_binned_cell_metadata_tmp <- data.frame()

    for (path in unique(processed_binned_cell_metadata[[path_colname]])) {
      # Get path
      binned.path.frame <- processed_binned_cell_metadata[processed_binned_cell_metadata[[path_colname]] == path, , drop = FALSE]

      # Generate refernce
      ref_time <- c(1:nrow(binned.path.frame))

      # Align with the binned Pseudotime
      binned.pseudotime <- binned.path.frame[[bin_pseudotime_colname]]

      # Check
      if (!all(ref_time == binned.pseudotime)) {
        binned.path.frame[[bin_pseudotime_colname]] <- ref_time
        tmp.bin.name <- paste0(path, "_bin_", ref_time)
        binned.path.frame[[bin_colname]] <- tmp.bin.name
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
      group_by(!!sym(path_colname)) %>%
      summarize(max_time = max(!!sym(bin_pseudotime_colname))) %>%
      ungroup()

    # Find the minimum of the max 'scmp_binned_pseudotime' across all 'Path'
    min_max_time <- min(pB.frame_tmp$max_time)

    # Identify rows to be removed
    rows_to_remove <- pB.frame %>%
      filter(!!sym(bin_pseudotime_colname) > min_max_time)

    # Print a message about the rows that will be removed
    if (verbose) {
      if (nrow(rows_to_remove) > 0) {
        message(paste("Dropped trailing bin", rows_to_remove[[bin_pseudotime_colname]], "from", rows_to_remove[[path_colname]]))
      }
    }

    # Filter out rows where 'scmp_binned_pseudotime' is greater than 'min_max_time'
    pB.frame <- pB.frame %>%
      filter(!!sym(bin_pseudotime_colname) <= min_max_time)

    processed_binned_cell_metadata <- pB.frame
  }

  compressed.sparse <- SingleCellExperiment(assays = list(bulk.counts = as(matrix(NA, nrow = 0, ncol = nrow(processed_binned_cell_metadata)), "dgCMatrix")))
  compressed.sparse@colData <- DataFrame(processed_binned_cell_metadata)
  scmpObject@dense <- compressed.sparse

  # Get Counts
  scmpObject <- pb_counts(
    scmpObject = scmpObject,
    bin_members_colname = bin_members_colname,
    bin_colname = bin_colname,
    assay_name = assay_name,
    cluster_count_by = cluster_count_by
  )

  # Update Slots
  scmpObject@param@pseudotime_colname <- pseudotime_colname
  scmpObject@param@path_colname <- path_colname
  scmpObject@param@bin_method <- bin_method
  scmpObject@param@bin_pseudotime_colname <- bin_pseudotime_colname
  scmpObject@param@bin_colname <- bin_colname
  scmpObject@param@bin_members_colname <- bin_members_colname
  scmpObject@param@bin_size_colname <- bin_size_colname
  return(scmpObject)
}
