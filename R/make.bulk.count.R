#' @title Create Pseduo-bulk Counts
#'
#' @description
#' `make.pseudobulk.counts()` creates a dataframe of pseudo bulk counts from single cell counts. It does this by either taking the mean or sum of counts across clusters in each bin, depending on the specified method.
#'
#' @param counts Raw count data from single cell RNA-seq experiment.
#' @param bin_members_colname Name of the column in the 'compressed_cell_metadata', 
#' storing information about the members of the bins. (Default is 'scmp_bin_members')
#' @param bin_colname Name of the column in the 'compressed_cell_metadata', 
#' storing information about the bin labels. (Default is 'scmp_bin')
#' @param pseudo_bulk_profile A data.frame generated using
#' \code{\link{make.pseudobulk.design}}.
#' @param cluster_count_by A character string (default = "sum"). The method to
#' use to aggregate counts within each cluster.
#'
#' @return
#' A data.frame. The data frame includes pseudo bulk counts with each row being a gene and each column being a bin.
#'
#' @details
#' The function operates by iterating over each row of the pseudo_bulk_profile. For each bin, it identifies the cells that belong to the bin and selects their counts from the counts data frame. It then calculates the mean or sum of these counts (depending on the specified method), and adds these to a new data frame of pseudo bulk counts. The result is a pseudo bulk counts data frame where each row is a gene and each column is a bin.
#'
#' @examples
#' \dontrun{
#' make.pseudobulk.counts(
#'   counts = sc_counts,
#'   pseudo_bulk_profile = pb_profile,
#'   cluster_count_by = "mean"
#' )
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' 
#' @export

make.pseudobulk.counts <- function(counts,
                                   bin_members_colname = "scmp_bin_members",
                                   bin_colname = "scmp_bin",
                                   pseudo_bulk_profile,
                                   cluster_count_by = "sum") {
  assert_that(bin_members_colname %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_members_colname, "' does not exist in pseudo-bulk-profile")
  )
  assert_that(bin_colname %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_colname, "' does not exist in pseudo-bulk-profile")
  )

  # Get the meta-information for pseudobulking
  meta.info <- pseudo_bulk_profile[, c(bin_members_colname, bin_colname)]

  # Determine the number of cores to use for parallel processing.
  # Here, I've used one less than the total number of cores available on the machine,
  # but you can adjust this based on your specific hardware.
  # num_cores <- detectCores() - 1

  # Run mclapply
  pb.counts <- lapply(1:nrow(meta.info), function(i) {
    # Get the bin.info
    bin <- meta.info[i, , drop = FALSE]

    # Split the row
    cell.vector <- c(str_split(bin[1], "\\|"))[[1]]

    # Get col cells
    col_indices <- which(colnames(counts) %in% cell.vector)

    # Subset the matrix using these indices
    bin_matrix <- as.matrix(counts[, col_indices, drop = F])

    # Get Pseudobulked-counts
    pb.vector <- switch(cluster_count_by,
      "mean" = as.matrix(round(rowMeans(bin_matrix))),
      "sum"  = as.matrix(rowSums(bin_matrix)),
      stop("Invalid cluster_count_by value. Please choose either 'mean' or 'sum'.")
    )

    # Return
    return(pb.vector)
  })

  # Convert the list output of mclapply to a matrix and set the row names
  pb.counts <- do.call(cbind, pb.counts)
  rownames(pb.counts) <- rownames(counts)
  colnames(pb.counts) <- meta.info[[bin_colname]]

  # return
  return(pb.counts)
}
