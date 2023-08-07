#' Create design matrix for 'scMaSigProClass' object
#'
#' This function creates a design matrix using the 'compress.sce' slot of a 'scMaSigProClass' object.
#' It generates an 'edesignClass' object which is then stored in the 'edesign' slot of the 'scMaSigProClass' object.
#'
#' @param scmpObj A 'scMaSigProClass' object.
#' @param degree Degree of the design matrix (default: 2).
#' @param time.col Name of the time column.
#' @param path.col Name of the path column.
#'
#' @return Returns the 'scmpObj' with an updated 'edesign' slot.
#' @export
#'
#' @examples
#' # Insert an example of how to use the function here.
#'
sc.make.design.matrix <- function(scmpObj,
                                  degree = 2,
                                  time.col,
                                  path.col) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Extract cell metadata
  comp.cell.metadata <- as.data.frame(colData(scmpObj@compress.sce))

  # Subset cell metadata
  com.cell.meta <- comp.cell.metadata[, colnames(comp.cell.metadata) %in% c(time.col, path.col)]

  # Get available paths
  avail.paths <- as.vector(unique(com.cell.meta[[path.col]]))

  # Add Dummy Variables
  for (i in avail.paths) {
    com.cell.meta[[i]] <- ifelse(com.cell.meta[[path.col]] %in% i, 1, 0)
  }

  # Drop path columns
  com.cell.meta <- com.cell.meta[, colnames(com.cell.meta) != path.col, drop = F]

  # Get colvec
  col.vec <- colnames(com.cell.meta)[colnames(com.cell.meta) != time.col]

  # Add Replicate Column
  com.cell.meta <- add_replicate_column(com.cell.meta, col.vec)

  # Order
  ord <- c(c(1, ncol(com.cell.meta)), c(2:c(ncol(com.cell.meta) - 1)))

  # Reorder columns
  com.cell.meta <- com.cell.meta[, ord]

  # Run Original MaSigPro make.matrix.design
  edesignList <- make.design.matrix(com.cell.meta,
    degree = degree,
    time.col = 1,
    repl.col = 2,
    group.cols = c(3:ncol(com.cell.meta))
  )

  # Create Object
  edesignObj <- new("edesignClass",
    dis = edesignList$dis,
    groups.vector = edesignList$groups.vector,
    edesign = edesignList$edesign
  )

  # Update Slot
  scmpObj@edesign <- edesignObj

  return(scmpObj)
}

add_replicate_column <- function(df, columns) {
  df <- df %>%
    mutate(Replicate = with(
      rle(paste0(.[, columns], collapse = "")),
      rep(seq_along(lengths), lengths)
    ))
  return(df)
}
