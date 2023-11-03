#' @title Pseudo-bulking with optimal number of pseudotime based bins
#'
#' @description
#' `squeeze()` discretizes a continuous time series column into bins
#' of equal size using entropy-based binning method. It automatically calculates
#' the optimal number of bins using one of the supported methods. The bin sizes
#' are also calculated and merged with the input cell_metadata.
#'
#' @param scmpObject object of Class scMaSigPro. See \code{\link{scMaSigProClass}}
#' for more details.
#' @param pseudotime_colname Name of the column in `cell.metadata` storing
#' Pseudotime values. Generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is "Pseudotime").
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path. Generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is `path_prefix`).
#' @param bin_pseudotime_colname Name of the column to store the computed Pseudotime
#' bins.
#' @param bin_colname Name of the bin column name
#' @param bin_size_colname Setting the name of the bin size column.
#' @param bin_method A character string specifying the method to use in order to
#' estimate the optimal number of bins. Available options: "Freedman.Diaconis",
#' "Sqrt", "Sturges", "Rice", "Doane", and "Scott.Normal". See \code{\link{estBinSize}}
#' for more details. (Default = "Sturges").
#' @param drop.fac A numeric value specifying the factor by which to decrease the
#' number of bins if the initial binning results in too many bins. (Default = 0.5).
#' @param verbose Print detailed output in the console. (Default is TRUE)
#' per path iteratively. Options: "universal", "individual. (Default = "universal").
#' @param homogenize_bins TRUE
#' @param additional_params Pass additional parameters as a named list. See Details.
#' @param assay_name Name of the Assay in the assay_name object from which retrieve the counts.
#' (Default = "counts").
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
#'
#' @examples
#' \dontrun{
#' squeeze(
#'   cell_metadata = data.frame, pseudotime_colname = "time",
#'   bin_method = "Sturges", drop.fac = 0.5, verbose = TRUE
#' )
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom parallel mclapply detectCores
#' @importFrom entropy discretize
#' @importFrom dplyr left_join join_by mutate select bind_rows group_by_at summarise rename_with
#' @importFrom magrittr %>%
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso \code{\link{estBinSize}}, \code{\link{discretize}}, \code{\link{create_range}}
#'
#' @export

squeeze <- function(scmpObject,
                    pseudotime_colname = scmpObject@addParams@pseudotime_colname,
                    path_colname = scmpObject@addParams@path_colname,
                    bin_method = "Sturges",
                    drop.fac = 0.5,
                    verbose = TRUE,
                    bin_colname = "scmp_bin",
                    bin_size_colname = "scmp_bin_size",
                    bin_pseudotime_colname = "scmp_binned_pseudotime",
                    homogenize_bins = TRUE,
                    assay_name = "counts",
                    additional_params = list(use_unique_time_points = FALSE)) {
  # Initiate Variable
  scmp_bin_lower_bound <- "scmp_l_bound"
  scmp_bin_upper_bound <- "scmp_u_bound"

  # Check Object Validity
  assert_that(is(scmpObject, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Extract cell metadata
  cell_metadata <- as.data.frame(colData(scmpObject@sce))

  # Drop Columns if exist
  cols_to_drop <- c(
    scmpObject@addParams@bin_size_colname,
    scmpObject@addParams@bin_pseudotime_colname,
    "scmp_u_bound", "scmp_l_bound"
  )
  cell_metadata <- cell_metadata[, !colnames(cell_metadata) %in% cols_to_drop, drop = FALSE]

  # Count slot
  assert_that(
    all(
      assay_name %in% names(scmpObject@sce@assays@data@listData)
    ),
    msg = paste0("'", assay_name, "' ", "doesn't exit in scmpObject.")
  )

  # Checks
  assert_that(pseudotime_colname %in% colnames(cell_metadata),
    msg = paste0("'", pseudotime_colname, "' does not exist in cell_metadata. Please review the 'pseudotime_colname' parameter.")
  )
  assert_that(path_colname %in% colnames(cell_metadata),
    msg = paste0("'", path_colname, "' does not exist in cell_metadata. Please review the 'path_colname' parameter.")
  )
  assert_that(drop.fac >= 0.3,
    msg = "Invalid value for 'drop.fac'. It should be between 0.3 and 1."
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
  cell_metadata$cell <- rownames(cell_metadata)

  # Get the avaible paths
  avail.paths <- as.vector(unique(cell_metadata[[path_colname]]))

  # Check for path
  assert_that(length(avail.paths) >= 2,
    msg = "Invalid number of paths detected. Please make sure that dataset has at least two paths"
  )

  if (verbose) {
    message(paste("Computing optimal bin-size with", bin_method, "method."))
    message(paste("Number of available path in the dataset:", length(avail.paths)))
    message(paste("Paths:", paste(avail.paths, collapse = ", ")))
    message(paste("Drop factor:", drop.fac))
  }

  # Apply transformations on data
  discrete.list <- lapply(avail.paths, function(path, design.frame = cell_metadata,
                                                drop_fac = drop.fac, path.col = path_colname,
                                                bin.size = bin_size_colname, bin = bin_colname,
                                                time.col = pseudotime_colname, method.bin = bin_method,
                                                bin.time.col = bin_pseudotime_colname,
                                                hom_bin = homogenize_bins,
                                                v = verbose, use.unique.time.points = additional_params$use_unique_time_points,
                                                lbound = scmp_bin_lower_bound, ubound = scmp_bin_upper_bound) {
    # Get the cells belonging to path
    path.frame <- design.frame[design.frame[[path.col]] == path, , drop = F]

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
        estBins <- estBinSize(
          time_vector = time_vector, nPoints = length_n,
          drop_fac = drop.fac, bin_method = method.bin
        )

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

    # Calculate Bin intervals with entropy
    bin_intervals <- as.data.frame(discretize(time_vector, numBins = estBins, r = range(time_vector)))
    
    # Client-Verbose
    if (verbose) {
      message(paste(
        "For", path, ",", length_n, "time points has been compressed to", nrow(bin_intervals), "bins"
      ))
    }

    # Clean the table before merge
    colnames(bin_intervals) <- c(bin, bin.size)
    
    # Create the bin table
    bin_table <- as.data.frame(
        t(
            as.data.frame(
                apply(
                    bin_intervals, 1, create_range,
                    bin_size_colname = bin.size,
                    bin_colname = bin, verbose = v
    ))))
    
    # Set column names
    colnames(bin_table) <- c(lbound, ubound, bin.size)
    
    # Get Mean and SD
    mean_value <- round(mean(bin_table[[bin.size]])/2)
    sd_value <- round(sd(bin_table[[bin.size]])/4)
    
    # Get thresholds
    max_allowed <- abs(mean_value + sd_value)
    min_allowed <- abs(mean_value - sd_value)
    
    print(paste("max_allowed:", max_allowed))
    print(paste("min_allowed:", min_allowed))
    
    # New Bin Map
    bin_table_uniform <- data.frame(matrix(NA, ncol = ncol(bin_table)))
    colnames(bin_table_uniform) <- colnames(bin_table)
    
    
    print("_-------------Raw----bin_table")
    
    print(bin_table)
    
    # Binning one-by one
    for (i in c(1:nrow(bin_table))){
        
        # histo_bin
        hist_bin <- c(unlist(bin_table[i,,drop = FALSE]))
        
        # Get size of the bins
        current_bin_size <- hist_bin[[bin.size]]
        
        print(paste("Current Size:", current_bin_size))
        
        # Check if the bin size is big or small
        if(current_bin_size > max_allowed){
            print("Greater than the max_allowed limit")
            
            potential_splits <- round(current_bin_size/max_allowed)
            print(paste("Number of potential bins to split", potential_splits))
            
            # get pseudotime
            pTime.per.interval <- time_vector[(time_vector >= hist_bin[[lbound]] & time_vector <= hist_bin[[ubound]])]
            
            # New Range
            new_range <- as.data.frame(entropy::discretize(pTime.per.interval, numBins = potential_splits, r = range(pTime.per.interval)))
            colnames(new_range) <- c(bin, bin.size)
            # Convert to table
            new_bin_table <- as.data.frame(
                t(as.data.frame(
                    apply(
                        new_range, 1, create_range,
                        bin_size_colname = bin.size,
                        bin_colname = bin, 
                        verbose = F
                    ))))
            # Set column names
            colnames(new_bin_table) <- c(lbound, ubound, bin.size)
            # Combine
            bin_table_uniform <- rbind(bin_table_uniform, new_bin_table)
        }else if(current_bin_size < min_allowed){
            print("Smaller than the max_allowed limit")
            
            # Get the prev and next bin
            prev_bin_size <- bin_table[i-1, ][[bin.size]]
            next_bin_size <- bin_table[i+1, ][[bin.size]]
            
            # Evaluate
            if(i <= 1 && prev_bin_size > next_bin_size & prev_bin_size >= min_allowed){
                print("Fuse from the previous bin")
            }else if(i > 1 & next_bin_size > prev_bin_size & next_bin_size >= min_allowed){
                print("Fuse from the next bin")
                
                print(hist_bin)
                next_bin_time_vector <- time_vector[(time_vector >= hist_bin[[lbound]])]
                next_bin_time_vector <- next_bin_time_vector[(next_bin_time_vector >= hist_bin[[ubound]])]
                
                print("------ordered_next_bin_time_vector-------")
                
                ordered_next_bin_time_vector <- sort(next_bin_time_vector)
                
                print(ordered_next_bin_time_vector)
                
                slice_first <- min_allowed + 1
                
                current_bin_new_elements <- ordered_next_bin_time_vector[c(1:slice_first)]
                
                print("------current_bin_new_elements-------")
                
                print(current_bin_new_elements)
                
                next_bin_new_elements <- ordered_next_bin_time_vector[c(slice_first:length(ordered_next_bin_time_vector))]
                
                
                print("------next_bin_new_elements-------")
                print(next_bin_new_elements)
                
                # Add to the dataframe
                new_range_current <- as.data.frame(entropy::discretize(current_bin_new_elements, numBins = 1, r = range(current_bin_new_elements)))
                new_range_next <- as.data.frame(entropy::discretize(next_bin_new_elements, numBins = 1, r = range(next_bin_new_elements)))
                colnames(new_range_current) <- c(bin, bin.size)
                colnames(new_range_next) <- c(bin, bin.size)
                # Convert to table
                new_bin_table_current <- as.data.frame(
                    t(as.data.frame(
                        apply(
                            new_range_current, 1, create_range,
                            bin_size_colname = bin.size,
                            bin_colname = bin, 
                            verbose = F
                        ))))
                # Set column names
                colnames(new_bin_table_current) <- c(lbound, ubound, bin.size)
                new_bin_table_next <- as.data.frame(
                    t(as.data.frame(
                        apply(
                            new_range_next, 1, create_range,
                            bin_size_colname = bin.size,
                            bin_colname = bin, 
                            verbose = F
                        ))))
                # Set column names
                colnames(new_bin_table_next) <- c(lbound, ubound, bin.size)
                
                print("------current_bin_new_elements-------")
                print(new_bin_table_current)
                
                
                print("------next_bin_new_elements-------")
                print(new_bin_table_next)
                
                
                bin_table_uniform <- rbind(bin_table_uniform, new_bin_table_current)
                
                bin_table[i+1,] <- new_bin_table_next
                
                print(bin_table_uniform)
                
                
                print("_-------------Replaced----bin_table")
                
                print(bin_table)
                }
            
        }else if(current_bin_size <= max_allowed & current_bin_size >= min_allowed){
            bin_table_uniform <- rbind(bin_table_uniform, hist_bin)
            print("Bin Size exist within the limits")
        }
    }
    bin_table_uniform <- bin_table_uniform[-1, ]
    
    if (path == "Path2"){
        print(bin_table)
        print(bin_table_uniform)
        stop()
    }
    
    if (verbose) {
      message(paste(
        "Finally, for", path, ",", length_n, "time points has been compressed to", nrow(bin_table), "bins"
      ))
    }

    # Create an empty data frame to store the results
    processed_cell_metadata <- data.frame()

    # Loop over each row of cell_metadata
    for (i in 1:nrow(path.frame)) {
      pseudotime_value <- path.frame[i, pseudotime_colname]

      # Filter bin_table based on pseudotime_value
      matching_bins <- bin_table[
        pseudotime_value >= bin_table[, lbound] &
          pseudotime_value <= bin_table[, ubound],
      ]

      # Combine the cell_metadata row with each matching bin
      combined_rows <- merge(path.frame[i, ], matching_bins, all.x = TRUE, by = character(0))

      # Bind combined_rows to results data frame
      processed_cell_metadata <- rbind(processed_cell_metadata, combined_rows)
    }
    # Convert result to a data frame
    processed_cell_metadata <- as.data.frame(processed_cell_metadata)

    # Set the 'cell' column as rownames
    rownames(processed_cell_metadata) <- processed_cell_metadata$cell

    return(processed_cell_metadata)
  })

  # Bind rows and convert to data frame, then drop 'cell' column
  processed_cell_metadata <- bind_rows(discrete.list) %>%
    as.data.frame()

  # Set the 'cell' column as rownames
  rownames(processed_cell_metadata) <- processed_cell_metadata$cell


  # Now, you can remove the 'cell' column
  processed_cell_metadata <- processed_cell_metadata %>% select(-"cell")

  ## Add Processed Cell Matadata back with slot update
  scmpObject@sce@colData <- DataFrame(processed_cell_metadata)

  # Update Slots
  scmpObject@addParams@pseudotime_colname <- pseudotime_colname
  scmpObject@addParams@path_colname <- path_colname
  scmpObject@addParams@bin_method <- bin_method
  scmpObject@addParams@bin_pseudotime_colname <- bin_pseudotime_colname
  scmpObject@addParams@bin_colname <- bin_colname
  scmpObject@addParams@bin_size_colname <- bin_size_colname
  return(scmpObject)
}
