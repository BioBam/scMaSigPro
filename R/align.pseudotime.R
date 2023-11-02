#' @importFrom scales rescale

align.pseudotime <- function(scmpObj, pseudotime_col, path_col, method = "rescale", verbose=TRUE){
    
    # Extract Cell metadata
    cell.metadata <- scmpObj@sce@colData %>% as.data.frame()
    
    # Extract Time and Path info with cell
    cell.metadata.sub <- cell.metadata[, c(pseudotime_col, path_col), drop=FALSE]
    cell.metadata.sub$cell <- rownames(cell.metadata.sub)
    
    # Get paths
    path.vec <- unique(cell.metadata.sub[[path_col]])
    
    # Get paths
    path1_time <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[1], pseudotime_col]
    names(path1_time) <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[1], "cell"]
    path2_time <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[2], pseudotime_col]
    names(path2_time) <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[2], "cell"]
    
    # Get list
    pTimeVectors <- select_longer_vector(vector1 = path1_time, vector1_label = path.vec[1],
                                         vector2 = path2_time, vector2_label = path.vec[2])
    
    if(length(pTimeVectors) == 4){
    # Update
    pTimeVectors$long_vec <- rescale(pTimeVectors$long_vec, to = c(min(pTimeVectors$short_vec), max(pTimeVectors$short_vec)),)
    
    short_tmp <- data.frame(time = pTimeVectors$short_vec,
                            cell = names(pTimeVectors$short_vec))
    colnames(short_tmp) <- c(paste(pseudotime_col, "rescaled", sep = "_"), "cell")
    long_tmp <- data.frame(time = pTimeVectors$long_vec,
                           cell = names(pTimeVectors$long_vec))
    colnames(long_tmp) <- c(paste(pseudotime_col, "rescaled", sep = "_"), "cell")
    new_time  <- rbind(short_tmp, long_tmp)
    
    # Merge
    cell.metadata <- merge(cell.metadata, new_time, by = "cell")
    rownames(cell.metadata) <- cell.metadata[["cell"]]
    
    # Drop cell
    cell.metadata <- cell.metadata[, colnames(cell.metadata) != "cell", drop = F]
    
    # Add
    scmpObj@sce@colData <- DataFrame(cell.metadata)
    
    # Update Pseudotime
    scmpObj@addParams@pseudotime_colname <- paste(pseudotime_col, "rescaled", sep = "_")
    
    return(scmpObj)
    }else{
        return(scmpObj)
    }
}
