#' @title Set up the polynomial design
#'
#' @description
#' `sc.set.poly()` builds the design of the maSigPro model (Conesa et al.
#' 2006): a 0/1 dummy for every path but the first, the reference, and the
#' products of each dummy with each power of the binned pseudotime up to
#' `poly_degree`.
#'
#' The reference is `reference` when it is given. Otherwise it is the first
#' path of the `Dense` slot: after [sc.squeeze()], the first level of a factor
#' path column, or else the first label in C-locale (byte) order; after
#' [sc.restruct()], the first path of `end_node_list`; with bins given to
#' [create_scmp()] or [cDense()], the path of the first bin. The other paths
#' keep the order of the `Dense` slot.
#'
#' @param poly_degree The degree of the polynomial in pseudotime.
#' @param reference The reference path, one of the path labels of the `Dense`
#'   slot, or `NULL` for the first path.
#'
#' @return The object with the design in the `Design` slot
#' (class [MatrixDesign-class]).
#'
#' @seealso [showPoly()], [predictors()] and [pathAssign()], which read the
#' design; [MatrixDesign-class].
#'
#' @references
#' Conesa A, Nueda MJ, Ferrer A, Talón M (2006). maSigPro: a method to
#' identify significantly differential expression profiles in time-course
#' microarray experiments. Bioinformatics 22(9), 1096-1102.
#' \doi{10.1093/bioinformatics/btl056}
#'
#' Srivastava P, Benegas Coll M, Götz S, Nueda MJ, Conesa A (2024).
#' scMaSigPro: differential expression analysis along single-cell
#' trajectories. Bioinformatics 40(7), btae443.
#' \doi{10.1093/bioinformatics/btae443}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#' @family workflow steps
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' scmp <- sc.set.poly(scmp.ob, poly_degree = 2)
#' showPoly(scmp)
#'
#' @export
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
    # Why: one 0/1 column per path; .poly_design() drops that of the
    # reference, so the terms of another path measure its difference
    # from the reference.
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

    # Create Object, without offset until a fit computes one
    designObj <- new("MatrixDesign",
        predictor_matrix = as.matrix(designList$dis),
        groups.vector = designList$groups.vector,
        assignment_matrix = as.matrix(designList$edesign)
    )

    # Update Slot
    scmpObj@Design <- designObj

    # Update poly degree
    scmpObj@Parameters@poly_degree <- as.integer(poly_degree)
    validObject(scmpObj)

    return(scmpObj)
}
