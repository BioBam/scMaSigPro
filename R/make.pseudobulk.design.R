#' @title Create Pseduo-bulk metadata
#'
#' @description
#' `make.bulk.design()` creates a pseudobulk profile data frame from the design file.
#' It groups cells into bins for each path from a given set of paths and
#' calculates the bin size.
#'
#' @param scmpObject object of Class scMaSigPro. See \code{\link{scMaSigProClass}}
#' for more details.
#' @param path_colname Name of the column in `cell.metadata` storing information
#' for Path. Generated using pkg{SingleCellExperiment}. (Default
#' is `path_prefix`).
#' @param bin_colname Name of the column in the 'compressed_cell_metadata',
#' storing information about the bin labels. (Default is 'scmp_bin').
#' @param bin_size_colname Name of the column in the 'compressed_cell_metadata'
#' storing information about the size of the bins. (Default is 'scmp_bin_size').
#' @param bin_members_colname Name of the column in the 'compressed_cell_metadata'
#' storing information about the members of the bins. (Default is 'scmp_bin_members').
#' @param bin_pseudotime_colname Name of the column in the 'compressed_cell_metadata'
#' storing information about the binned pseudotime. (Default is 'scmp_binned_pseudotime').
#' @param fill_gaps description
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return
#' A data.frame containing the pseudobulk profile. The data frame includes the following columns:
#' - 'binnedTime': The time bin.
#' - 'cluster.members': The cells grouped in the bin.
#' - 'bin': The bin identifier.
#' - 'path': The path identifier.
#' - 'bin.size': The size of the bin (number of cells in the bin).
#' - Columns for each path in paths.vector, with binary values indicating whether the row belongs to the path.
#'
#' @details
#' This function operates by iterating over the specified paths. For each path, it:
#' - Filters the design file to only include cells that belong to the path.
#' - Groups the cells into bins based on the binned time column.
#' - Calculates the size of each bin.
#' - Appends the resulting data frame to the pseudobulk profile.
#'
#' After processing all paths, the function adds binary columns for each path to the pseudobulk profile.
#'
#' @examples
#' \dontrun{
#' make.pseudobulk.design(
#'   compressed_cell_metadata = df,
#'   path_colname = "path1", binned_pseudotime_column = "binnedTime"
#' )
#' }
#'
#' @seealso \code{\link{calc_bin_size}}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom rlang :=
#' @importFrom dplyr group_by filter row_number ungroup summarise summarize
#'
#' @export
make.pseudobulk.design <- function(scmpObject,
                                   path_colname = scmpObject@addParams@path_colname,
                                   bin_colname = scmpObject@addParams@bin_colname,
                                   bin_size_colname = scmpObject@addParams@bin_size_colname,
                                   bin_members_colname = "scmp_bin_members",
                                   bin_pseudotime_colname = scmpObject@addParams@bin_pseudotime_colname,
                                   verbose = TRUE,
                                   fill_gaps = FALSE) {
  # Check Object Validity
  assert_that(is(scmpObject, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Extract cell metadata
  compressed_cell_metadata <- as.data.frame(colData(scmpObject@sce))
  assert_that(bin_pseudotime_colname %in% colnames(compressed_cell_metadata),
    msg = paste0("'", bin_pseudotime_colname, "' does not exist in compressed_cell_metadata, please run entropy_discretize()")
  )
  assert_that(path_colname %in% colnames(compressed_cell_metadata),
    msg = paste0("'", path_colname, "' does not exist in compressed_cell_metadata. Please review 'path_colname' parameter.")
  )

  # Initate Variable
  scmp_bar <- "scmp_bar"

  # Get the avaible paths
  avail.paths <- as.vector(unique(compressed_cell_metadata[[path_colname]]))

  # Add helper-col
  compressed_cell_metadata[[scmp_bar]] <- rownames(compressed_cell_metadata)

  # Check for path
  assert_that(length(avail.paths) >= 2,
    msg = "Invalid number of paths detected. Please make sure that dataset has atleast two paths"
  )

  # Determine the number of cores
  # num_cores <- detectCores() - 1 # leave one core free for other tasks

  # Apply transformations on data
  # pB.list <- mclapply(avail.paths, function(path, design.frame = compressed_cell_metadata,
  pB.list <- lapply(avail.paths, function(path, design.frame = compressed_cell_metadata,
                                          binned.col = bin_pseudotime_colname, path.col = path_colname,
                                          v = verbose, fill.gaps = fill_gaps) {
    if (v) {
      message("Running for", path)
    }
    # Get the cells belonging to path
    path.frame <- design.frame[design.frame[[path.col]] == path, , drop = F]

    # Order along the temporal vector
    path.time.cell <- path.frame[order(path.frame[, binned.col]), c(binned.col, scmp_bar)]

    # Validation
    assert_that(nrow(path.time.cell) >= 2,
      msg = paste("Time points are already less/equal than/to two in", path)
    )

    # Group by time
    path.time.cell <- path.time.cell %>%
      group_by_at(binned.col) %>%
      summarise(!!sym(bin_members_colname) := paste0(scmp_bar, collapse = "|"))

    # Add Cluster Label
    tmp.bin.name <- paste0(path, "_bin_", seq(1, nrow(path.time.cell)))
    path.time.cell[[bin_colname]] <- tmp.bin.name

    # Set the Path Information
    path.time.cell[[path.col]] <- path

    # Add Cluster Size
    tmp.bin.sise <- apply(path.time.cell, 1, calc_bin_size, clus_mem_col = bin_members_colname)
    path.time.cell[[bin_size_colname]] <- tmp.bin.sise

    if (v) {
      message(paste(tmp.bin.name, "has", tmp.bin.sise, "cells."))
    }

    # Claculate bin_range
    bin_range <- range(path.time.cell[[bin_size_colname]], na.rm = TRUE)

    # throw warning
    warningCondition(
      diff(bin_range) >= 100,
      paste("Differences among bin sizes are greater than 100 units for", path)
    )

    if (fill.gaps) {
      # if fill groups are true
      path.time.cell[[binned.col]] <- rownames(path.time.cell)
    }

    # Return frame
    return(path.time.cell)
    # }, mc.cores = num_cores)
  })

  # Bind rows
  pB.frame <- bind_rows(pB.list) %>% as.data.frame()

  # Add rownames
  rownames(pB.frame) <- pB.frame[[bin_colname]]

  # Remove extra column
  # pB.frame <- pB.frame %>% select(-"scmp_bar")

  if (fill_gaps) {
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
  }

  # Add rownames
  rownames(pB.frame) <- pB.frame[[bin_colname]]

  # Covert columns to mumeric
  pB.frame[[bin_pseudotime_colname]] <- as.numeric(pB.frame[[bin_pseudotime_colname]])

  ## Add Processed Cell Matadata back with slot update
  compressed.sce <- SingleCellExperiment(assays = list(bulk.counts = as(matrix(NA, nrow = 0, ncol = nrow(pB.frame)), "dgCMatrix")))
  compressed.sce@colData <- DataFrame(pB.frame)
  scmpObject@compress.sce <- compressed.sce

  ## Slot Update
  scmpObject@addParams@path_colname <- path_colname
  scmpObject@addParams@bin_colname <- bin_colname
  scmpObject@addParams@bin_size_colname <- bin_size_colname
  scmpObject@addParams@bin_members_colname <- bin_members_colname
  scmpObject@addParams@bin_pseudotime_colname <- bin_pseudotime_colname

  # Pathway infor
  return(scmpObject)
}
