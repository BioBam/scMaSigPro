#' make.pseudobulk.counts
#'
#' @description
#' This function creates a dataframe of pseudo bulk counts from single cell counts. It does this by either taking the mean or sum of counts across clusters in each bin, depending on the specified method.
#'
#' @param counts A data.frame. The count data from single cell RNA-seq experiment.
#' @param cluster_member_col A character string (default = "cluster.members"). The name of the column in the pseudo_bulk_profile that contains the cluster members.
#' @param bin_col A character string (default = "bin"). The name of the column in the pseudo_bulk_profile that contains the bin identifiers.
#' @param pseudo_bulk_profile A data.frame. The pseudo bulk profile created by make_pseudobulk_design function.
#' @param cluster.count.by A character string (default = "sum"). The method to use to aggregate counts within each cluster. Possible values are "mean" or "sum".
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
#'   cluster.count.by = "mean"
#' )
#' }
#'
#' @export

make.pseudobulk.counts <- function(counts, cluster_member_col = "cluster.members",
                                   bin_col = "bin", pseudo_bulk_profile,
                                   cluster.count.by = "sum") {
  assert_that(cluster_member_col %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", cluster_member_col, "' does not exist in pseudo-bulk-profile")
  )
  assert_that(bin_col %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_col, "' does not exist in pseudo-bulk-profile")
  )

  # Get the meta-information for pseudobulking
  meta.info <- pseudo_bulk_profile[, c(cluster_member_col, bin_col)]

  # Define the function to be used in mclapply
  compress <- function(bin, count.table = counts, cluster.count.by) {
    cell.vector <- c(str_split(bin[1], "\\|"))[[1]]

    # Find column indices that match cell.vector
    col_indices <- which(colnames(count.table) %in% cell.vector)

    # Subset the matrix using these indices
    bin_matrix <- count.table[, col_indices]

    # Get Pseudobulked-counts
    pb.vector <- switch(cluster.count.by,
      "mean" = as.matrix(rowMeans(bin_matrix)),
      "sum"  = as.matrix(rowSums(bin_matrix)),
      stop("Invalid cluster.count.by value. Please choose either 'mean' or 'sum'.")
    )
    return(pb.vector)
  }

  # Determine the number of cores to use for parallel processing.
  # Here, I've used one less than the total number of cores available on the machine,
  # but you can adjust this based on your specific hardware.
  num_cores <- detectCores() - 1

  # Run mclapply
  pb.counts <- mclapply(1:nrow(meta.info), function(i) compress(meta.info[i, , drop = FALSE], count.table = counts, cluster.count.by = cluster.count.by), mc.cores = num_cores)

  # Convert the list output of mclapply to a matrix and set the row names
  pb.counts <- do.call(cbind, pb.counts)
  rownames(pb.counts) <- rownames(counts)
  colnames(pb.counts) <- meta.info[[bin_col]]

  # return
  return(pb.counts)
}
