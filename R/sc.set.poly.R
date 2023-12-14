#' Create predictors and set polynomial.  Adaption of maSigPro::make.design.matrix()
#'
#' @param scmpObject A 'scmp' object.
#' @param poly_degree Degree of the design matrix (default: 2).
#' @param bin_pseudotime_colname Name of the time column.
#' @param path_colname Name of the path column.
#'
#' @return Returns the 'scmpObject' with an updated 'design' slot.
#'
#' @importFrom maSigPro make.design.matrix
#' @export
#'
sc.set.poly <- function(scmpObject,
                        poly_degree = 2,
                        bin_pseudotime_colname = scmpObject@param@bin_pseudotime_colname,
                        path_colname = scmpObject@param@path_colname) {
  # Check Object Validity
  assert_that(is(scmpObject, "scmp"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Extract cell metadata
  comp.cell.metadata <- as.data.frame(scmpObject@dense@colData)

  # pseudotime_colname
  assert_that((bin_pseudotime_colname %in% colnames(comp.cell.metadata)),
    msg = paste0("'", bin_pseudotime_colname, "' ", "doesn't exit in cell.metadata.")
  )
  # path_colname
  assert_that((path_colname %in% colnames(comp.cell.metadata)),
    msg = paste0("'", path_colname, "' ", "doesn't exit in cell.metadata.")
  )

  # Subset cell metadata
  com.cell.meta <- comp.cell.metadata[, colnames(comp.cell.metadata) %in% c(bin_pseudotime_colname, path_colname)]

  # Get available paths
  avail.paths <- as.vector(unique(com.cell.meta[[path_colname]]))

  # Add Dummy Variables
  for (i in avail.paths) {
    com.cell.meta[[i]] <- ifelse(com.cell.meta[[path_colname]] %in% i, 1, 0)
  }

  # Drop path columns
  com.cell.meta <- com.cell.meta[, colnames(com.cell.meta) != path_colname, drop = FALSE]

  # Get colvec
  col.vec <- colnames(com.cell.meta)[colnames(com.cell.meta) != bin_pseudotime_colname]

  # Add Replicate Column
  # com.cell.meta <- com.cell.meta %>%
  #   mutate(Replicate = data.table::rleid(Reduce(paste, com.cell.meta))) %>%
  #   as.data.frame()
  com.cell.meta$Replicate <- with(com.cell.meta, {
    # Create a concatenated string of all columns
    combined <- apply(com.cell.meta, 1, paste, collapse = "-")
    # Use rle (run length encoding) to find runs of identical values
    rle_ids <- with(rle(combined), rep(seq_along(lengths), lengths))
    return(rle_ids)
  })

  # Order
  ord <- c(c(1, ncol(com.cell.meta)), c(2:c(ncol(com.cell.meta) - 1)))

  # Reorder columns
  com.cell.meta <- as.matrix(com.cell.meta[, ord])

  # Run Original MaSigPro make.matrix.design
  designList <- make.design.matrix(com.cell.meta,
    degree = poly_degree,
    time.col = 1,
    repl.col = 2,
    group.cols = c(3:ncol(com.cell.meta))
  )

  # Create Object
  designObj <- new("designClass",
    predictor = as.matrix(designList$dis),
    groups.vector = designList$groups.vector,
    alloc = as.matrix(designList$edesign)
  )

  # Update Slot
  scmpObject@design <- designObj

  # Update poly degree
  scmpObject@param@poly_degree <- as.integer(poly_degree)

  return(scmpObject)
}
