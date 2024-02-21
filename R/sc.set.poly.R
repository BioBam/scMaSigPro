#' @title Set up polynomial models and create Predictor Matrix
#'
#' @description
#' Set up polynomial models and create Predictor Matrix that will contain the
#' independent variables. It is a wrapper around `maSigPro::make.design.matrix`.
#'
#' @importFrom maSigPro make.design.matrix
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param poly_degree Degree of the polynomial.
#' @param bin_ptime_col A character string representing the column name
#' for binned Pseudotime values in 'Dense' data.
#' @param path_col A character string representing the column name for branching
#' path assignment in 'Sparse' or 'Dense' slot.
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with updated `Design`
#' slot.
#'
#' @seealso \code{\link{MatrixDesign}} Class.
#'
#' @references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles
#' in Time-Course Microarray Experiments. Bioinformatics 22, 1096-1102}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#' @export
#'
sc.set.poly <- function(scmpObj,
                        poly_degree = 2,
                        bin_ptime_col = scmpObj@Parameters@bin_ptime_col,
                        path_col = scmpObj@Parameters@path_col) {
  # Check Object Validity
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Extract cell metadata
  comp.cell.metadata <- as.data.frame(scmpObj@Dense@colData)

  # pseudotime_colname
  assert_that((bin_ptime_col %in% colnames(comp.cell.metadata)),
    msg = paste0("'", bin_ptime_col, "' ", "doesn't exit in cell.metadata.")
  )
  # path_col
  assert_that((path_col %in% colnames(comp.cell.metadata)),
    msg = paste0("'", path_col, "' ", "doesn't exit in cell.metadata.")
  )

  # Subset cell metadata
  com.cell.meta <- comp.cell.metadata[, colnames(comp.cell.metadata) %in% c(bin_ptime_col, path_col)]

  # Get available paths
  avail.paths <- as.vector(unique(com.cell.meta[[path_col]]))

  # Add Dummy Variables
  for (i in avail.paths) {
    com.cell.meta[[i]] <- ifelse(com.cell.meta[[path_col]] %in% i, 1, 0)
  }

  # Drop path columns
  com.cell.meta <- com.cell.meta[, colnames(com.cell.meta) != path_col, drop = FALSE]

  # Get colvec
  col.vec <- colnames(com.cell.meta)[colnames(com.cell.meta) != bin_ptime_col]

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

  # Create offset
  offset_vector <- numeric(rep(nrow(designList$edesign)))
  names(offset_vector) <- rownames(designList$edesign)

  # Create Object
  designObj <- new("MatrixDesign",
    predictor_matrix = as.matrix(designList$dis),
    groups.vector = designList$groups.vector,
    assignment_matrix = as.matrix(designList$edesign),
    offset = offset_vector
  )

  # Update Slot
  scmpObj@Design <- designObj

  # Update poly degree
  scmpObj@Parameters@poly_degree <- as.integer(poly_degree)

  return(scmpObj)
}
