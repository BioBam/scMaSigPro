#' @title Restructure the bins of paths that share a root
#'
#' @description
#' `sc.restruct()` rebuilds the bins so that every end path starts at the
#' root: the bins of each link path are added to the end paths it links, the
#' bins of the root to every end path, the binned pseudotime of the later
#' bins is shifted to follow them, and [pb_counts()] recomputes the
#' pseudo-bulk counts.
#'
#' @param end_node_list A list of the labels of the end paths.
#' @param root_node The label of the root path.
#' @param link_node_list A list of the labels of the link paths, the segments
#'   that some end paths share; a label holds the labels of the paths it
#'   links, separated by `link_sep`.
#' @param link_sep The separator in the labels of the link paths, matched as
#'   fixed text when the label contains it and as a regular expression
#'   otherwise, so `"|"` and `"\\|"` both work.
#'
#' @return The object with the bins of every end path, from the root on,
#' and their pseudo-bulk counts in the `Dense` slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @seealso [assign_paths()], which labels the root and link paths.
#'
#' @family workflow steps
#' @inheritParams scmp-params
#'
#' @examples
#' data("multi.lin.sce", package = "scMaSigPro")
#' scmp <- as_scmp(multi.lin.sce,
#'     additional_params = list(
#'         labels_exist = TRUE, exist_ptime_col = "Monocle3_Pseudotime",
#'         exist_path_col = "group"
#'     ),
#'     verbose = FALSE
#' )
#' scmp <- sc.squeeze(scmp)
#' # Three end paths, one link path shared by two of them, and the root
#' scmp <- sc.restruct(scmp,
#'     end_node_list = list("path_Y_15", "path_Y_18", "path_Y_52"),
#'     root_node = "root", link_node_list = list("path_Y_18|path_Y_15"),
#'     link_sep = "|", verbose = FALSE
#' )
#' table(cDense(scmp)$Path)
#'
#' @export
sc.restruct <- function(scmpObj,
                        end_node_list,
                        root_node,
                        link_node_list,
                        assay_name = "counts",
                        aggregate = "sum",
                        link_sep = "_links_",
                        verbose = TRUE) {
    # Check Object Validity
    .check_scmp(scmpObj)
    .check_strings(assay_name = assay_name)
    aggregate <- .check_choice(aggregate, c("sum", "mean"),
        msg = "Invalid 'aggregate' value. Please choose either 'mean' or 'sum'."
    )

    # Check if values are binned
    .check(nrow(as.data.frame(colData(scmpObj@Dense))) >= 1,
        msg = "No binning information found. Please run 'sc.squeeze()', first."
    )

    # Count slot
    .check(
        all(
            assay_name %in% assayNames(scmpObj@Sparse)
        ),
        msg = paste0("'", assay_name, "' ", "doesn't exist in scmpObj.")
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
        # Extract link group
        link_value <- link_node_list[[i]]
        names(link_value) <- link_value

        # Split on link_sep as fixed text when the label contains it, so
        # that "|" works, else as a regular expression, as "\\|" needs
        fixed_sep <- isTRUE(grepl(link_sep[1L], link_value, fixed = TRUE))
        link_vec <- unlist(strsplit(link_value, link_sep, fixed = fixed_sep))

        # Traverse for each of the link
        for (j in link_vec) {
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
    linked <- bin_info[[scmpObj@Parameters@path_col]] %in%
        bin_link_tmp[[scmpObj@Parameters@path_col]]
    binfo_tmp <- bin_info[!linked, ]
    binfo_tmp[["scmp_restruct"]] <- "NA"
    bin_link_tmp <- rbind(bin_link_tmp, binfo_tmp)

    # Run end-point wise operation
    for (i in end_node_list) {
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
            message("Linking root bins to '", i, "'")
        }
    }

    new_bin_info <- bin_root_tmp

    # Bin metadata only; pb_counts() adds the assay
    scmpObj@Dense <- SingleCellExperiment::SingleCellExperiment(
        colData = S4Vectors::DataFrame(new_bin_info)
    )

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
