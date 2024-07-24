#' @title Restructure the binned data.
#'
#' @description
#' `sc.restruct()` Add Description
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' path assignment in 'Sparse' or 'Dense' data.
#' @param end_node_list A list of end nodes in of the branch.
#' @param root_node A character string specifying the root node.
#' @param link_node_list A list of links between two nodes.
#' @param assay_name Name of the Assay in sparse data from which the counts are
#' used. (Default = "counts").
#' @param aggregate A character string specifying the method to aggregate counts
#' within each cluster. Available options are 'mean' or 'sum'. (Default = "sum").
#' @param verbose  Print detailed output in the console. (Default is TRUE)
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with updated `Dense`
#' slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso \code{\link{estBinSize}}, \code{\link{discretize}},
#' \code{\link{create_range}}
#'
#' @export
sc.restruct <- function(scmpObj,
                        end_node_list,
                        root_node,
                        link_node_list,
                        assay_name = "counts",
                        aggregate = "sum",
                        verbose = TRUE) {
    # scmpObj <- multi_scmp_ob
    # root_node <- "root"
    # end_node_list <- list("path_Y_15", "path_Y_18", "path_Y_52")
    # link_node_list <- list("path_Y_18|path_Y_15")
    # assay_name <- "counts"
    # aggregate <- "sum"
    # verbose = TRUE
    
    # Check Object Validity
    assertthat::assert_that(is(scmpObj, "ScMaSigPro"),
                            msg = "Please provide object of class 'scMaSigPro'."
    )
    
    # Check if values are binned
    assertthat::assert_that(nrow(as.data.frame(colData(scmpObj@Dense))) >= 1,
                            msg = "No binning information found. Please run 'sc.squeeze()', first."
    )
    
    # Count slot
    assertthat::assert_that(
        all(
            assay_name %in% names(scmpObj@Sparse@assays@data@listData)
        ),
        msg = paste0("'", assay_name, "' ", "doesn't exit in scmpObj.")
    )
    
    # Extract bin information
    bin_info <- cDense(scmpObj)
    
    # Set list names with values
    names(link_node_list) <- unlist(link_node_list)
    names(end_node_list) <- unlist(end_node_list)
    
    # Invoke Empty dfs
    bin_link_tmp <- data.frame()
    bin_root_tmp <- data.frame()
    
    # Add Links to paths
    for (i in link_node_list) {
        ## remove
        # i <- link_node_list[[1]]
        
        # Extract link group
        link_value <- link_node_list[[i]]
        names(link_value) <- link_value
        
        # Split
        link_vec <- unlist(stringr::str_split(link_value, "\\|"))
        
        # Traverse for each of the link
        for (j in link_vec) {
            ## remove
            # j <- link_vec[1]
            
            # Extract end_bin_info
            link_end_bin_info <- bin_info[bin_info[[scmpObj@Parameters@path_col]] == j, ]
            
            # Extract link_bin_info
            link_bin_info <- bin_info[bin_info[[scmpObj@Parameters@path_col]] == i, ]
            
            # Update group of link_bin_info
            link_bin_info[[scmpObj@Parameters@path_col]] <- j
            
            # Update rowname and bin_col
            new_label <- paste(link_bin_info[[scmpObj@Parameters@path_col]], "bin", link_bin_info[[scmpObj@Parameters@bin_ptime_col]], sep = "_")
            rownames(link_bin_info) <- new_label
            link_bin_info[[scmpObj@Parameters@bin_col]] <- new_label
            
            # Calculate offset
            link_offset <- max(link_bin_info[[scmpObj@Parameters@bin_ptime_col]])
            
            # Add new column for scmp_restruct
            link_bin_info[["scmp_restruct"]] <- i
            
            # Add offset to end_bin_info
            link_end_bin_info[[scmpObj@Parameters@bin_ptime_col]] <- link_end_bin_info[[scmpObj@Parameters@bin_ptime_col]] + link_offset
            
            # Update label and bin_col
            new_label <- paste(link_end_bin_info[[scmpObj@Parameters@path_col]], "bin", link_end_bin_info[[scmpObj@Parameters@bin_ptime_col]], sep = "_")
            link_end_bin_info[[scmpObj@Parameters@bin_col]] <- new_label
            rownames(link_end_bin_info) <- new_label
            
            # Add restructure column
            link_end_bin_info[["scmp_restruct"]] <- j
            
            # Combine
            tmp <- rbind(link_bin_info, link_end_bin_info)
            
            # Rbind
            bin_link_tmp <- rbind(bin_link_tmp, tmp)
            
            
            if (verbose) {
                message("Linking bins from '", i, "' to path '", j, "'")
            }
        }
    }
    
    # Create end
    names(end_node_list) <- unlist(end_node_list)
    
    # Extract Root bins
    root_node_info <- bin_info[bin_info[[scmpObj@Parameters@path_col]] %in% root_node, ]
    root_node_info[["scmp_restruct"]] <- rep("root", nrow(root_node_info))
    root_offset <- max(root_node_info[[scmpObj@Parameters@bin_ptime_col]])
    
    # Get rows
    binfo_tmp <- bin_info[!(bin_info[[scmpObj@Parameters@path_col]] %in% bin_link_tmp$group), ]
    binfo_tmp[["scmp_restruct"]] <- "NA"
    bin_link_tmp <- rbind(bin_link_tmp, binfo_tmp)
    
    # Run end-point wise operation
    for (i in end_node_list) {
        ## remove
        # leaf <- end_node_list[[1]]
        
        # Extract end_bin_info for leaf
        leaf <- end_node_list[[i]]
        
        # Extract leaf bin info
        leaf_bin_info <- bin_link_tmp[bin_link_tmp[[scmpObj@Parameters@path_col]] == leaf, ]
        
        # Add Offset
        leaf_bin_info[[scmpObj@Parameters@bin_ptime_col]] <- leaf_bin_info[[scmpObj@Parameters@bin_ptime_col]] + root_offset
        
        # Calculate new label
        new_label <- paste(leaf_bin_info[[scmpObj@Parameters@path_col]], "bin", leaf_bin_info[[scmpObj@Parameters@bin_ptime_col]], sep = "_")
        
        # Update
        rownames(leaf_bin_info) <- new_label
        leaf_bin_info[[scmpObj@Parameters@bin_col]] <- new_label
        leaf_bin_info[["scmp_restruct"]] <- leaf
        
        # Update roor bin info
        root_node_info_tmp <- root_node_info
        root_node_info_tmp[[scmpObj@Parameters@path_col]] <- leaf
        new_label <- paste(root_node_info_tmp[[scmpObj@Parameters@path_col]], "bin", root_node_info_tmp[[scmpObj@Parameters@bin_ptime_col]], sep = "_")
        
        # Update
        rownames(root_node_info_tmp) <- new_label
        root_node_info_tmp[[scmpObj@Parameters@bin_col]] <- new_label
        
        # Add
        tmp <- rbind(root_node_info_tmp, leaf_bin_info)
        bin_root_tmp <- rbind(bin_root_tmp, tmp)
        
        if (verbose) {
            message("Linking bins root bins to '", i)
        }
    }
    
    new_bin_info <- bin_root_tmp
    
    compressed.sparse <- SingleCellExperiment::SingleCellExperiment(assays = list(
        bulk.counts = as(
            matrix(NA, nrow = 0, ncol = nrow(new_bin_info)),
            "dgCMatrix"
        )
    ))
    
    compressed.sparse@colData <- S4Vectors::DataFrame(new_bin_info)
    scmpObj@Dense <- compressed.sparse
    
    # Get Counts
    scmpObj <- pb_counts(
        scmpObj = scmpObj,
        bin_mem_col = scmpObj@Parameters@bin_mem_col,
        bin_col = scmpObj@Parameters@bin_col,
        assay_name = assay_name,
        cluster_count_by = aggregate
    )
    return(scmpObj)
}
