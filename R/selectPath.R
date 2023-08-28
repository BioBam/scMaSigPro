#' Select Paths for scMaSigPro
#'
#' This function selects paths to be used in scMaSigPro. It extracts
#' a sub-object based on the specified paths.
#'
#' @param obj An object of class `scMaSigProClass`. This object will be checked
#'   to ensure it's the right type.
#' @param sel.path A character vector indicating the paths to be selected.
#' @param pathCol The column name in `colData(sceObj)` that indicates the path.
#'   Defaults to "Path".
#'
#' @return A `SingleCellExperiment` object subsetted based on the specified paths.
#'
#' @examples
#' # Assuming you have an example object of class `scMaSigProClass`:
#' # selected_obj <- selectPath(example_obj, sel.path = c("Path1", "Path3"))
#'
#' @export
selectPath <- function(obj, sel.path,
                       pathCol = "Path"){
    
    # Check
    assert_that(is(obj)[1] == "scMaSigProClass",
                msg = "Please supply object from scMaSigPro Class")
    
    # Extract the sce class
    sceObj <- obj@sce
    
    # Select Cells
    sceObj_sub <- sceObj[, colData(sceObj)[[pathCol]] %in% sel.path]
    
    # Add the Object Back
    obj@sce <- sceObj_sub
    
    # return
    return(obj)
}
