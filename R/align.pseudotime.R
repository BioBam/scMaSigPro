align.pseudotime <- function(scmpObj, pseudotime_col, path_col, method = "dtw", verbose=TRUE){
    
    # Extract Cell metadata
    cell.metadata <- scmpObj@sce@colData %>% as.data.frame()
    
    # Extract Time and Path info with cell
    cell.metadata.sub <- cell.metadata[, c(pseudotime_col, path_col), drop=FALSE]
    cell.metadata.sub$cell <- rownames(cell.metadata.sub)
    
    # Extract unique paths
    paths <- unique(cell.metadata.sub[, path_col])
    if (length(paths) != 2) {
        stop("The function currently supports exactly 2 paths.")
    }
    
    pseudotime_path1 <- cell.metadata.sub[cell.metadata.sub[, path_col] == paths[1], pseudotime_col]
    pseudotime_path2 <- cell.metadata.sub[cell.metadata.sub[, path_col] == paths[2], pseudotime_col]
    
    if (method == "simple_rescale") {
        # Use the simple rescale method
        rescaled_values <- rescale_pseudotime(pseudotime_path1, pseudotime_path2, verbose)
        cell.metadata.sub[cell.metadata.sub[, path_col] == paths[1], pseudotime_col] <- rescaled_values$path1
        cell.metadata.sub[cell.metadata.sub[, path_col] == paths[2], pseudotime_col] <- rescaled_values$path2
    } else {
        # Compute the range difference for each path
        diff_range_path1 <- diff(range(pseudotime_path1))
        diff_range_path2 <- diff(range(pseudotime_path2))
        
        # Rescale the pseudotimes of the path with the larger range
        if (diff_range_path1 > diff_range_path2) {
            if (verbose) cat("Rescaling pseudotime of", paths[1], "to match the range of", paths[2], "\n")
            # Rescale pseudotime of path1 to fit into the range of path2
            scaling_factor <- diff_range_path2 / diff_range_path1
            rescaled_pseudotime_path1 <- min(pseudotime_path2) + scaling_factor * (pseudotime_path1 - min(pseudotime_path1))
            cell.metadata.sub[cell.metadata.sub[, path_col] == paths[1], pseudotime_col] <- rescaled_pseudotime_path1
        } else {
            if (verbose) cat("Rescaling pseudotime of", paths[2], "to match the range of", paths[1], "\n")
            # Rescale pseudotime of path2 to fit into the range of path1
            scaling_factor <- diff_range_path1 / diff_range_path2
            rescaled_pseudotime_path2 <- min(pseudotime_path1) + scaling_factor * (pseudotime_path2 - min(pseudotime_path2))
            cell.metadata.sub[cell.metadata.sub[, path_col] == paths[2], pseudotime_col] <- rescaled_pseudotime_path2
        }
    }
    
    # Add to the metadata
    new_pseudotime <- paste(pseudotime_col, "dtw", sep = "_")
    colnames(cell.metadata.sub) <- c(new_pseudotime, path_col, "cell")
    
    # leftjoin over two columns
    cell.metadata$cell <- rownames(cell.metadata)
    cell.metadata.dtw <- merge(cell.metadata.sub, cell.metadata, by.x = c("cell", path_col), by.y = c("cell", path_col), all.x = TRUE)
    rownames(cell.metadata.dtw) <- cell.metadata.dtw$cell
    cell.metadata.dtw <- cell.metadata.dtw[, !colnames(cell.metadata.dtw) %in% "cell"]
    
    scmpObj@sce@colData <- DataFrame(cell.metadata.dtw)
    scmpObj@addParams@pseudotime_colname <- new_pseudotime
    
    # Return the modified data
    return(scmpObj)
}

rescale_pseudotime <- function(path1, path2, verbose = TRUE) {
    # Step 1: Calculate the range for both paths
    range_path1 <- range(path1)
    range_path2 <- range(path2)
    
    # Step 2: Check which path has a smaller range
    if (verbose) {
        cat("Range of path1:", range_path1, "\n")
        cat("Range of path2:", range_path2, "\n")
    }
    
    if (diff(range_path1) < diff(range_path2)) {
        # Rescale path2 to range of path1
        if (verbose) cat("Rescaling path2 to range of path1.\n")
        rescaled_path2 <- rescale_to_range(path2, range_path1)
        return(list(path1 = path1, path2 = rescaled_path2))
    } else {
        # Rescale path1 to range of path2
        if (verbose) cat("Rescaling path1 to range of path2.\n")
        rescaled_path1 <- rescale_to_range(path1, range_path2)
        return(list(path1 = rescaled_path1, path2 = path2))
    }
}

rescale_to_range <- function(data, new_range) {
    old_range <- range(data)
    scale_factor <- diff(new_range) / diff(old_range)
    rescaled_data <- (data - old_range[1]) * scale_factor + new_range[1]
    return(rescaled_data)
}
