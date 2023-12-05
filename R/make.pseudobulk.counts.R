#' @title Create Pseduo-bulk Counts
#'
#' @description
#' `make.pseudobulk.counts()` creates a dataframe of pseudo bulk counts from single
#' cell counts. It does this by either taking the mean or sum of counts across clusters
#' in each bin, depending on the specified method.
#'
#' @param scmpObject object of Class scMaSigPro. See \code{\link{scMaSigProClass}}
#' for more details.
#' @param bin_members_colname Column name in the 'compressed_cell_metadata'
#' storing information about the members of the bins. (Default is 'scmp_bin_members').
#' @param bin_colname Column name in the 'compressed_cell_metadata'
#' storing information about the bin labels. (Default is 'scmp_bin').
#' @param cluster_count_by A character string specifying the method to use to
#' aggregate counts within each cluster. Available options are 'mean' or 'sum'. (Default = "sum").
#' @param assay_name Name of the Assay in the assay_name object from which retrieve the counts.
#' (Default = "counts").
#'
#' @return
#' A data.frame. The data frame includes pseudo bulk counts with each row being
#' a gene and each column being a bin.
#'
#' @details
#' The function operates by iterating over each row of the pseudo_bulk_profile.
#' For each bin, it identifies the cells that belong to the bin and selects their
#' counts from the counts data frame. It then calculates the mean or sum of these
#' counts (depending on the specified method), and adds these to a new data frame of pseudo bulk counts.
#' The result is a pseudo bulk counts data frame where each row is a gene and each column is a bin.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @keywords internal

make.pseudobulk.counts <- function(scmpObject,
                                   bin_members_colname = scmpObject@param@bin_members_colname,
                                   bin_colname = scmpObject@param@bin_colname,
                                   assay_name = "counts",
                                   cluster_count_by = "sum") {
  # Check Object Validity
  assert_that(is(scmpObject, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'."
  )

  # Count slot
  assert_that(
    all(
      assay_name %in% names(scmpObject@sce@assays@data@listData)
    ),
    msg = paste0("'", assay_name, "' ", "doesn't exit in scmpObject.")
  )

  # Get assay
  counts <- scmpObject@sce@assays@data@listData[[assay_name]]

  # Get Pseudobulk Profile
  pseudo_bulk_profile <- as.data.frame(colData(scmpObject@compress.sce))

  assert_that(bin_members_colname %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_members_colname, "' does not exist in level.meta.data")
  )
  assert_that(bin_colname %in% colnames(pseudo_bulk_profile),
    msg = paste0("'", bin_colname, "' does not exist in level.meta.data")
  )

  # Get the meta-information for pseudobulking
  meta.info <- pseudo_bulk_profile[, c(bin_members_colname, bin_colname)]

  # Run mclapply
  pb.counts <- lapply(1:nrow(meta.info), function(i) {
    # Get the bin.info
    bin <- meta.info[i, , drop = FALSE]

    # Split the row
    cell.vector <- c(str_split(bin[1], "\\|"))[[1]]

    # Get col cells
    col_indices <- which(colnames(counts) %in% cell.vector)

    # Subset the matrix using these indices
    bin_matrix <- as.matrix(counts[, col_indices, drop = FALSE])

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

  # Return the counts
  scmpObject@compress.sce@assays@data@listData$bulk.counts <- as(pb.counts, "dgCMatrix")

  # return
  return(scmpObject)
}
