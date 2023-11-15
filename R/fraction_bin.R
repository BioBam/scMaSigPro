#' Binning and Plotting Function for scmpObj Data
#'
#' This function takes a single-cell compression object (scmpObj) and performs 
#' data processing, binning, and plotting. It extracts necessary information 
#' from the compressed and cell data within the object, processes this data for 
#' visualization, and creates a ggplot object representing the distribution 
#' of cells across different bins and paths.
#'
#' @param scmpObj An object representing compressed single-cell data. It is 
#'                expected to contain specific nested data structures and 
#'                parameters used for processing.
#'
#' @return A ggplot object representing the binned data in a bar plot format.
#'         The plot shows the distribution of cells across different bins, 
#'         time points, and paths.
#'
#' @examples
#' # Assuming 'myScmpObj' is a pre-existing single-cell compression object
#' plot <- sc.fraction.bin(myScmpObj)
#' print(plot)
#'
#' @importFrom ggplot2 ggplot geom_bar scale_x_continuous labs theme_minimal theme element_text
#' @importFrom dplyr mutate group_by summarise
#' @importFrom magrittr %>%
#' @importFrom stringr str_split_1
#' @export
sc.fraction.bin <- function(scmpObj,
                            pseudotime_colname = scmpObj@addParams@pseudotime_colname,
                            path_colname = scmpObj@addParams@path_colname,
                            annotation_col = scmpObj@addParams@annotation_col) {
    # Extract necessary columns from compressed and cell data
    compressedData <- as.data.frame(scmpObj@compress.sce@colData)
    cellData <- as.data.frame(scmpObj@sce@colData)
    
    # Add cell identifier to cell data
    cellDataSubset <- cellData[, c(pseudotime_colname, 
                                   path_colname, 
                                   annotation_col), drop = FALSE]
    cellDataSubset[["cell"]] <- rownames(cellDataSubset)

    # Process compressed data
    barPlotData <- lapply(1:nrow(compressedData), function(i) {
        rowVector <- compressedData[i, , drop = FALSE]
        cells <- str_split_1(rowVector[[scmpObj@addParams@bin_members_colname]], "\\|")
        data.frame(
            bin = rep(rowVector[[scmpObj@addParams@bin_colname]], length(cells)),
            row_time = rep(rowVector[[scmpObj@addParams@bin_pseudotime_colname]], length(cells)),
            cell = cells,
            path = rep(rowVector[[scmpObj@addParams@path_colname]], length(cells))
        )
    }) %>% do.call("rbind", .)
    
    # Merge with cell data to include additional information
    mergedData <- merge(cellDataSubset, barPlotData, by = "cell")
    mergedData$row_time <- as.numeric(mergedData$row_time)
    
    # Prepare data for plotting
    mergedData <- mergedData %>%
        mutate(row_time = factor(row_time, levels = unique(row_time)),
               Path = factor(Path, levels = c("Path1", "Path2")),
               cell_type = factor(!!sym(annotation_col), levels = unique(!!sym(annotation_col))))
    
    # Create summary data
    dfSummary <- mergedData %>%
        group_by(row_time, Path, cell_type) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        mutate(interaction = as.numeric(row_time) + 
                   (as.numeric(Path) - 1) * (0.9 / length(unique(Path))))
    
    
    dfSummary$row_time <- as.numeric(dfSummary$row_time)
    #return(head(dfSummary))
    
    
    
    # Generate and return plot
    dodgeWidth <- 0.9 / length(unique(dfSummary$Path))
    dfSummary$row_time <- factor(dfSummary$row_time, levels = sort(unique(dfSummary$row_time)))
    p <- ggplot(dfSummary, aes(x = interaction, y = Count, fill = cell_type)) +
        geom_bar(stat = "identity", position = "stack", width = dodgeWidth, alpha  =1) +
        scale_x_continuous("Row Time and Path", 
                           breaks = seq_along(unique(dfSummary$row_time)), 
                           labels = levels(dfSummary$row_time)) +
        labs(y = "Count", fill = "Cell Type") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    return(p)
}
