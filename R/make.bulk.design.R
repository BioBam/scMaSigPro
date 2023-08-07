#' make.pseudobulk.design
#'
#' @description
#' This function creates a pseudobulk profile data frame from the design file. It groups cells into bins for each path from a given set of paths and calculates the bin size.
#'
#' @param design.file A data.frame. The design file containing the binned data.
#' @param pathCol A character string. The name of the column in the design file that contains the path data.
#' @param binnedCol A character string (default = "binnedTime"). The name of the column in the design file that contains the binned data.
#'
#' @return
#' A data.frame containing the pseudobulk profile. The data frame includes the following columns:
#' - 'binnedTime': The time bin
#' - 'cluster.members': The cells that fall into the bin
#' - 'bin': The bin identifier
#' - 'path': The path identifier
#' - 'bin.size': The size of the bin (number of cells in the bin)
#' - Columns for each path in paths.vector, with binary values indicating whether the row belongs to the path
#'
#' @details
#' This function operates by iterating over the specified paths. For each path, it:
#' - Filters the design file to only include cells that belong to the path
#' - Groups the cells into bins based on the binned time column
#' - Calculates the size of each bin
#' - Appends the resulting data frame to the pseudobulk profile
#'
#' After processing all paths, the function adds binary columns for each path to the pseudobulk profile.
#'
#' @examples
#' \dontrun{
#' make.pseudobulk.design(design.file = df, pathCol = "path1", binnedCol = "binnedTime")
#' }
#'
#' @seealso \code{\link{calc_bin_size}}
#'
#' @export

make.pseudobulk.design <- function(design.file, pathCol,
                                   binnedCol = "binnedTime") {
  assert_that(binnedCol %in% colnames(design.file),
    msg = paste0("'", binnedCol, "' does not exist in design.file")
  )
  assert_that(pathCol %in% colnames(design.file),
    msg = paste0("'", pathCol, "' does not exist in design.file")
  )

  # Get the avaible paths
  avail.paths <- as.vector(unique(design.file[[pathCol]]))

  # Add helper-col
  design.file$scmp_bar <- rownames(design.file)

  # Check for path
  assert_that(length(avail.paths) >= 2,
    msg = "Invalid number of paths detected. Please make sure that dataset has atleast two paths"
  )

  # Determine the number of cores
  num_cores <- detectCores() - 1 # leave one core free for other tasks

  # Apply transformations on data
  pB.list <- mclapply(avail.paths, function(path, design.frame = design.file,
                                            # pB.list <- lapply(avail.paths, function(path, design.frame = design.file,
                                            binned.col = binnedCol, path.col = pathCol) {
    # Get the cells belonging to path
    path.frame <- design.frame[design.frame[[path.col]] == path, , drop = F]

    # Order along the temporal vector
    path.time.cell <- path.frame[order(path.frame[, binned.col]), c(binned.col, "scmp_bar")]

    # Validation
    assert_that(nrow(path.time.cell) >= 2,
      msg = paste("Time points are already less/equal than/to two in", path)
    )

    # Group by time
    path.time.cell <- path.time.cell %>%
      group_by_at(binned.col) %>%
      summarise(cluster.members = paste0(scmp_bar, collapse = "|"))

    # Add Cluster Label
    path.time.cell$bin <- paste0(path, "_bin_", seq(1, nrow(path.time.cell)))

    # Set the Path Information
    path.time.cell$path <- path

    # Add Cluster Size
    path.time.cell$bin.size <- apply(path.time.cell, 1, calc_bin_size)

    # Claculate bin_range
    bin_range <- range(path.time.cell$bin.size, na.rm = TRUE)

    # throw warning
    warningCondition(
      diff(bin_range) >= 100,
      paste("Differences among bin sizes are greater than 100 units for", path)
    )

    # Return frame
    return(path.time.cell)
  }, mc.cores = num_cores)
  # })

  # Bind rows
  pB.frame <- bind_rows(pB.list) %>% as.data.frame()

  # Add rownames
  rownames(pB.frame) <- pB.frame$bin

  # Remove extra column
  # pB.frame <- pB.frame %>% select(-"scmp_bar")

  # Pathway infor
  return(pB.frame)
}
