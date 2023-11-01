#' @title Find optimal number of pseudotime bins
#'
#' @description
#' `entropy_discretize()` discretizes a continuous time series column into bins
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
#' @param binning A character string. When set to "individual", the bins are calculated
#' per path iteratively. Options: "universal", "individual. (Default = "universal").
#' @param additional_params Pass additional parameters as a named list. See Details.
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
#' entropy_discretize(
#'   cell_metadata = data.frame, pseudotime_colname = "time",
#'   bin_method = "Sturges", drop.fac = 0.5, verbose = TRUE
#' )
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom parallel mclapply detectCores
#' @importFrom entropy discretize
#' @importFrom dplyr left_join join_by mutate select bind_rows group_by_at summarise rename_with closest
#' @importFrom magrittr %>%
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso \code{\link{estBinSize}}, \code{\link{discretize}}, \code{\link{create_range}}
#'
#' @export

entropy_discretize <- function(scmpObject,
                               pseudotime_colname = scmpObject@addParams@pseudotime_colname,
                               path_colname = scmpObject@addParams@path_colname,
                               bin_method = "Sturges",
                               drop.fac = 0.5,
                               verbose = TRUE,
                               binning = "universal",
                               bin_colname = "scmp_bin",
                               bin_size_colname = "scmp_bin_size",
                               bin_pseudotime_colname = "scmp_binned_pseudotime",
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
  assert_that(all(binning %in% c("universal", "individual")),
    msg = "Allowed options for binning are 'universal' and 'individual'"
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
    message(paste("Initiating binning by:", binning, "method."))
  }

  # Switch
  result <- switch(binning,
    universal = {
      # Extract the time information as a vector
      time_vector <- cell_metadata[, pseudotime_colname]
      length_n <- length(time_vector)

      if (additional_params$use_unique_time_points) {
        time_vector <- unique(time_vector)
        length_n <- length(time_vector)
        if (verbose) {
          message(paste("Using only unique points in the time series"))
        }
      }

      if (verbose) {
        message(paste("Number of pseudotime points detected", length_n))
        message(paste("Range of Pseudotime points", paste(range(time_vector), collapse = "-")))
      }

      # Calculate Optimal Number of Bins
      tryCatch(
        expr = {
          estBins <- estBinSize(
            time_vector = time_vector, nPoints = length_n,
            drop_fac = drop.fac, bin_method = bin_method
          )
          if (verbose) {
            message(paste("Sucessfully estimated optimal number of bin size."))
          }
        },
        error = function(e) {
          message(paste("Error message: ", e$message))
          stop("Unable to estimate bin size")
        }
      )

      if (verbose) {
        message(paste("Estimated bin size is", estBins))
      }

      # Calculate Bin intervals with entropy
      bin_intervals <- as.data.frame(discretize(time_vector, numBins = estBins, r = range(time_vector)))

      # Clean the table before merge
      colnames(bin_intervals) <- c(bin_colname, bin_size_colname)
      bin_intervals[[bin_pseudotime_colname]] <- rownames(bin_intervals)

      # Create the bin table
      bin_table <- as.data.frame(t(as.data.frame(apply(bin_intervals, 1, create_range,
        bin_pseudotime_colname = bin_pseudotime_colname,
        bin_size_colname = bin_size_colname, bin_colname = bin_colname, verbose = verbose
      ))))
      colnames(bin_table) <- c(scmp_bin_lower_bound, scmp_bin_upper_bound, bin_size_colname, bin_pseudotime_colname)

      if (verbose) {
        message(paste("Estimating bin intervals"))
      }

      # Combine Tables
      processed_cell_metadata <- as.data.frame(
        left_join(cell_metadata, bin_table,
          by = join_by(
            closest(!!pseudotime_colname >= !!scmp_bin_lower_bound),
            closest(!!pseudotime_colname <= !!scmp_bin_upper_bound)
          )
        )
      )

      # Set the 'cell' column as rownames
      rownames(processed_cell_metadata) <- processed_cell_metadata$cell
    },
    individual = {
      # Apply transformations on data
      discrete.list <- lapply(avail.paths, function(path, design.frame = cell_metadata,
                                                    drop_fac = drop.fac, path.col = path_colname,
                                                    bin.size = bin_size_colname, bin = bin_colname,
                                                    time.col = pseudotime_colname, method.bin = bin_method,
                                                    bin.time.col = bin_pseudotime_colname,
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
        bin_intervals[[bin.time.col]] <- rownames(bin_intervals)

        # Create the bin table
        bin_table <- as.data.frame(t(as.data.frame(apply(bin_intervals, 1, create_range,
          bin_pseudotime_colname = bin.time.col,
          bin_size_colname = bin.size, bin_colname = bin, verbose = v
        ))))
        colnames(bin_table) <- c(lbound, ubound, bin.size, bin.time.col)

        # Combine Tables
        processed_cell_metadata <- as.data.frame(
          left_join(path.frame, bin_table,
            by = join_by(
              closest(!!time.col >= !!lbound),
              closest(!!time.col <= !!ubound)
            )
          )
        )

        return(processed_cell_metadata)
      })

      # Bind rows and convert to data frame, then drop 'cell' column
      processed_cell_metadata <- bind_rows(discrete.list) %>%
        as.data.frame()

      # Set the 'cell' column as rownames
      rownames(processed_cell_metadata) <- processed_cell_metadata$cell
    },
    stop("Invalid option")
  )

  # Now, you can remove the 'cell' column
  processed_cell_metadata <- processed_cell_metadata %>% select(-"cell")

  ## Add Processed Cell Matadata back with slot update
  scmpObject@sce@colData <- DataFrame(processed_cell_metadata)

  # Update Slots
  scmpObject@addParams@pseudotime_colname <- pseudotime_colname
  scmpObject@addParams@path_colname <- path_colname
  scmpObject@addParams@bin_method <- bin_method
  scmpObject@addParams@binning <- binning
  scmpObject@addParams@bin_pseudotime_colname <- bin_pseudotime_colname
  scmpObject@addParams@bin_colname <- bin_colname
  scmpObject@addParams@bin_size_colname <- bin_size_colname
  return(scmpObject)
}
