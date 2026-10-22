#' @title Set up polynomial models and create Predictor Matrix
#'
#' @description
#' Set up polynomial models and create Predictor Matrix that will contain the
#' independent variables: the design of the maSigPro model (Conesa et al.
#' 2006), with a dummy for every path but the first, the reference, and its
#' products with each power of the binned pseudotime up to `poly_degree`.
#'
#' The reference is `reference` when it is given. Otherwise it is the first
#' path of the `Dense` slot: after [sc.squeeze()], the first level of the path
#' column when it is a factor and else the first path label in C-locale
#' (byte) order, whatever the locale of the session; after [sc.restruct()],
#' the first path of `end_node_list`; with bins given to [create_scmp()] or
#' \code{\link{cDense<-}}, the path of the first bin. The other paths keep the
#' order of the `Dense` slot.
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param poly_degree Degree of the polynomial.
#' @param bin_ptime_col A character string representing the column name
#' for binned Pseudotime values in 'Dense' data.
#' @param path_col A character string representing the column name for branching
#' path assignment in 'Sparse' or 'Dense' slot.
#' @param reference The reference path, one of the path labels of the `Dense`
#' slot, or `NULL` (default) for the first path of the `Dense` slot. See
#' Description.
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
                        path_col = scmpObj@Parameters@path_col,
                        reference = NULL) {
    # Check Object Validity
    .check_scmp(scmpObj)
    .check_strings(bin_ptime_col = bin_ptime_col, path_col = path_col)
    if (!is.null(reference)) {
        .check_strings(reference = reference)
    }

    # Extract cell metadata
    comp.cell.metadata <- as.data.frame(colData(scmpObj@Dense))

    # pseudotime_colname
    .check((bin_ptime_col %in% colnames(comp.cell.metadata)),
        msg = paste0(
            "'", bin_ptime_col, "' ", "doesn't exist in cell.metadata."
        )
    )
    # path_col
    .check((path_col %in% colnames(comp.cell.metadata)),
        msg = paste0("'", path_col, "' ", "doesn't exist in cell.metadata.")
    )

    # Subset cell metadata
    com.cell.meta <- comp.cell.metadata[, colnames(comp.cell.metadata) %in% c(bin_ptime_col, path_col)]

    # Get available paths: the reference first, then the order of the bins
    avail.paths <- as.vector(unique(com.cell.meta[[path_col]]))
    if (!is.null(reference)) {
        .check(reference %in% avail.paths,
            msg = paste0(
                "'reference' must be one of the paths: ",
                paste(avail.paths, collapse = ", "), "."
            )
        )
        avail.paths <- c(reference, setdiff(avail.paths, reference))
    }

    # Add Dummy Variables
    for (i in avail.paths) {
        com.cell.meta[[i]] <- ifelse(com.cell.meta[[path_col]] %in% i, 1, 0)
    }

    # Drop path columns
    com.cell.meta <- com.cell.meta[, colnames(com.cell.meta) != path_col, drop = FALSE]

    # Get colvec
    col.vec <- colnames(com.cell.meta)[colnames(com.cell.meta) != bin_ptime_col]

    # Add Replicate Column
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

    # Build the design: time in column 1, the paths from column 3
    designList <- .poly_design(com.cell.meta, degree = poly_degree)

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
    validObject(scmpObj)

    return(scmpObj)
}
