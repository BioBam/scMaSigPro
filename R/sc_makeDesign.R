#' Create Object for Sc MaSigPro
#'
#' @param sce.object Object of the \link[SingleCellExperiment]{SingleCellExperiment} class.
#' @param time.col Column name having pseudotime information.
#' @param path.col Column name having path information.
#' @param poly.order Max degree of the polynomial
#'

# Function to Create ScMaSigPro Object
sc_makeDesign <- function(sce.object, path.col = "Group",
                                  time.col = "Step", poly.order = 2) {
  # Extract the cell Metadata
  cell.meta <- as.data.frame(colData(sce.object))

  # Extract the levels
  path.vec.levels <- levels(as.factor(cell.meta[[path.col]]))
  groups <- path.vec.levels

  # Extract the values
  path.vec.values <- as.character(cell.meta[[path.col]])

  # Binarize and Create a design file
  encoded.frame <- as.matrix(sapply(path.vec.levels, function(X, path.values = path.vec.values) {
    encoded.vec <- ifelse(path.values == X, 1, 0)
    return(encoded.vec)
  }))

  # Add time Information
  time.values <- as.numeric(cell.meta[[time.col]])

  # Reverse the Vector, select base path and update path
  base.path <- path.vec.levels[1]
  path.vec.levels <- path.vec.levels[-1]

  # Creation of Vs_Covariate
  vs_covar <- sapply(path.vec.levels, function(X, base_path = base.path,
                                               encoded_frame = encoded.frame) {
    vs.covar <- paste(X, base_path, sep = "_vs_")
    vs.covar.value <- list(as.numeric(encoded_frame[, X, drop = TRUE]))
    names(vs.covar.value) <- vs.covar
    return(vs.covar.value)
  }, USE.NAMES = F)

  # Join the vectors
  vs_covar <- do.call(cbind, vs_covar)

  # Create Interaction with Time
  int_covar <- lapply(c(1:poly.order), function(X, encoded_frame = encoded.frame, time_col = time.col,
                                                time_values = time.values, path_vec_levels = path.vec.levels) {
    time.poly <- as.numeric(time.values**X)

    time.ploy.int <- sapply(path_vec_levels, function(Y, encoded_frame_tmp = encoded_frame, time_poly = time.poly) {
      path.dummys <- encoded_frame_tmp[, Y, drop = T]
      path.dummys <- as.numeric(path.dummys * time_poly)
      return(path.dummys)
    }, USE.NAMES = T)

    # Renaming for the list
    time.poly <- as.matrix(time.poly)
    colnames(time.poly) <- paste(time_col, X, sep = "")
    colnames(time.ploy.int) <- paste(colnames(time.poly), colnames(time.ploy.int), sep = "_x_")
    int.covar <- cbind(time.poly, time.ploy.int)
    return(int.covar)
  })

  # Combine the list
  int_covar <- do.call(cbind, int_covar)

  # Combine to make the covariate Matrix
  covariates <- cbind(vs_covar, int_covar)

  # Create Group Vector
  group_vector.tmp <- c(colnames(vs_covar), base.path)
  group_vector <- rep(group_vector.tmp, ncol(covariates))
  group_vector <- group_vector[c(1:ncol(covariates))]

  parameters.obj <- new("parameterClass",
    path.col = path.col,
    time.col = time.col,
    poly.order = as.integer(poly.order)
  )

  covar.obj <- new("covariateClass",
    factor.design = encoded.frame,
    time.series = time.values,
    covariate = covariates,
    covariate.vector = colnames(covariates),
    group.vector = group_vector,
    groups = c(groups)
  )

  # Add back as a new slot in the S4
  scMaSigPro.object <- new("scMaSigProClass", sce.object,
    covariate = covar.obj,
    parameters = parameters.obj
  )

  # Return object
  return(scMaSigPro.object)
}
