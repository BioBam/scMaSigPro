#' Select Paths for scMaSigPro
#'
#' This function selects paths to be used in scMaSigPro. It extracts
#' a sub-object based on the specified paths.
#'
#' @param obj An object of class `scMaSigProClass`. This object will be checked
#'   to ensure it's the right type.
#' @param sel.path A character vector indicating the paths to be selected.
#' @param balance_paths Logical indicating if paths should be balanced. Default is TRUE.
#' @param pathCol The column name in `colData(sceObj)` that indicates the path.
#'   Defaults to "Path".
#' @param plot_paths Logical indicating if paths should be plotted. Default is TRUE.
#' @param pTimeCol The column name indicating the pseudotime. Default is "Pseudotime".
#' @param verbose Logical indicating if verbose messages should be displayed. Default is TRUE.
#'
#' @return A `SingleCellExperiment` object subsetted based on the specified paths.
#'
#' @examples
#' # Assuming you have an example object of class `scMaSigProClass`:
#' # selected_obj <- selectPath(example_obj, sel.path = c("Path1", "Path3"))
#'
#' @export
selectPath <- function(obj, sel.path, balance_paths = T,
                       pathCol = "Path", plot_paths = T,
                       pTimeCol = "Pseudotime", verbose = T) {
    # Check
    assert_that(is(obj)[1] == "scMaSigProClass",
                msg = "Please supply object from scMaSigPro Class"
    )
    
    # Extract the sce class
    sceObj <- obj@sce
    
    # Extract Cell MetaData
    cell.meta.raw <- as.data.frame(colData(sceObj))
    
    # Extract Paths from the metadata
    cell.meta.sub <- cell.meta.raw[cell.meta.raw[[pathCol]] %in% sel.path, ]
    
    # Plot
    if (plot_paths) {
        before.plt <- ggplot(cell.meta.sub, aes(x = .data[[pathCol]], y = .data[[pTimeCol]])) +
            geom_boxplot(color = "#f58a53") + # This creates the boxplots for each category
            stat_summary(fun = median, geom = "line", aes(group = 1), color = "#e84258") +
            stat_summary(fun = median, geom = "point", color = "#159287") +
            labs(
                title = "A. Before Balance",
                x = "Available Paths", y = "Inferred Pseudotime"
            ) +
            theme_minimal(base_size = 15)
    }
    
    # Balance
    if (balance_paths) {
        cell.meta.sub <- slice_by_smallest_range(
            data = cell.meta.sub, path_col = pathCol, pseudotime_col = pTimeCol
        )
        
        if (plot_paths) {
            after.plt <- ggplot(cell.meta.sub, aes(x = .data[[pathCol]], y = .data[[pTimeCol]])) +
                geom_boxplot(color = "#e84258") + # This creates the boxplots for each category
                stat_summary(fun = median, geom = "line", aes(group = 1), color = "#f58a53") +
                stat_summary(fun = median, geom = "point", color = "#159287") +
                labs(
                    title = "B. After Balance",
                    x = "Available Paths", y = "Inferred Pseudotime"
                ) +
                theme_minimal(base_size = 15)
        }
    }
    
    # Plot
    if (plot_paths) {
        comb.plt <- ggarrange(before.plt, after.plt)
        print(comb.plt)
    }
    
    # Select Cells
    sceObj_sub <- sceObj[, rownames(colData(sceObj)) %in% rownames(cell.meta.sub)]
    
    # Add the Object Back
    obj@sce <- sceObj_sub
    
    # return
    return(obj)
}
