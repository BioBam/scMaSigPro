#' Slice Data By The Smallest Pseudotime Range
#'
#' This function trims data from paths based on the path with the smallest pseudotime range.
#' Specifically, it determines which path has the smallest range of pseudotimes
#' and then slices all other paths to match this range.
#'
#' @param data A dataframe containing single cell data with columns specified by `path_col` and `pseudotime_col`.
#' @param path_col A string specifying the column name in `data` that contains path identifiers.
#' @param pseudotime_col A string specifying the column name in `data` that contains pseudotime values.
#' @param verbose A logical value indicating whether to print informational messages. Default is `TRUE`.
#'
#' @return A dataframe where paths have been trimmed to match the pseudotime range of the path with the smallest range.
#' Cells outside this range for each path are removed.
#'
slice_by_smallest_range <- function(data, path_col, pseudotime_col, verbose = T) {
    
    # Detect unique paths
    paths <- unique(data[[path_col]])
    
    # Store original ranges
    original_ranges <- lapply(paths, function(path) {
        return(
            range(data[data[[path_col]] == path, pseudotime_col], na.rm = TRUE)
            )
    })
    
    # Name the list 
    names(original_ranges) <- paths
    
    # Get the span of the ranges
    range_spans <- sapply(original_ranges, diff)
    
    # Order the ranges 
    sorted_paths <- names(range_spans[order(range_spans)])
    
    # First path is the smallest
    smallest_path <- sorted_paths[1]
    
    if(!is.na(smallest_path)){
    
    # Get the range of the smallest path
    smallest_range <- original_ranges[[smallest_path]]
    
    # Set rownames to columns
    data$row_id <- rownames(data)
    
    # Filter
    filtered_data <- data %>%
        group_by(across(all_of(path_col))) %>%
        filter((!!sym(pseudotime_col) >= smallest_range[1]) & (!!sym(pseudotime_col) <= smallest_range[2])) %>%
        as.data.frame()
    
    # Optionally, set "row_id" as row names again
    rownames(filtered_data) <- filtered_data$row_id
    filtered_data$row_id <- NULL
    
    if (verbose) {
        removed_count <- nrow(data) - nrow(filtered_data)
        message(paste(removed_count, "cells were removed to match the pseudotime range of", smallest_path))
    }
    
    return(filtered_data)
    }else{
        message("Nothing to remove, paths correspond")
        return(data)
    }
    
}

