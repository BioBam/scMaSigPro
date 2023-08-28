#' Convert cell_data_set or SingleCellExperiment to scMaSigProClass
#'
#' This function converts a cell_data_set object from Monocle or a SingleCellExperiment object to an instance of
#' the scMaSigProClass. Currently, only Monocle's cell_data_set and SingleCellExperiment are supported.
#'
#' @param object An S4 object of class 'cell_data_set' or 'SingleCellExperiment'.
#' @param from Character string specifying the class of 'object'. Can be either "cell_data_set" or "sce".
#'
#' @return An instance of the 'scMaSigProClass'.
#' @examples
#' \dontrun{
#' scmpObj <- as_scmp(object, from = "sce")
#' }
#'
#' @importFrom monocle3 cell_data_set
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom assertthat assert_that
#' @export
as_scmp <- function(object, from = "cell_data_set", verbose = T) {
    
    # Check Conversion Type
    assert_that(from %in% c("cell_data_set", "sce"),
                msg = ("Currently, accepted formats are S4 cell_data_set and sce")
    )
    
    # Validate S4
    assert_that(
        all(isS4(object) & all(is(object, "cell_data_set") | is(object, "SingleCellExperiment"))),
        msg = "Please provide object from one of the class 'Monocle3/cell_data_set', 'SingleCellExperiment/SCE'"
    )
    
    # FLow control
    if (is(object)[1] == "SingleCellExperiment") {
        if(verbose){message("SCE/SingleCellExperiment Detected")}
        # Create Object
        scmpObj <- new("scMaSigProClass",
                       sce = object,
                       compress.sce = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
        )
        # Return Object
        return(scmpObj)
    } else if (is(object)[1] == "cell_data_set") {
        if(verbose){message("Monocle3/cell_data_set Detected")}
        
        # Annotate the monocel3 Object
        annotated_cds <- extract_monocle3_components(cds, reduction_method = "umap")
        
        # Convert to sce abd then to scMaSigPro Class
        scmpObj <- new("scMaSigProClass",
                       sce = annotated_cds,
                       compress.sce = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
        )
        return(scmpObj)
    }
}
