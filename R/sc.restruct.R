#' @title Restructure the bins of paths that share a root
#'
#' @description
#' `sc.restruct()` rebuilds the bins so that every end path starts at the
#' root: the bins of each link path are added to the end paths it links, the
#' bins of the root to every end path, the binned pseudotime of the later
#' bins is shifted to follow them, and [pb_counts()] recomputes the
#' pseudo-bulk counts.
#'
#' @param end_node_list A list of the labels of the end paths. The bins of
#'   a path that is not the root, an end path or a link path are left out,
#'   with a message when `verbose` is `TRUE`.
#' @param root_node The label of the root path. Each path is named once, in
#'   `root_node`, `end_node_list` or `link_node_list`.
#' @param link_node_list A list of the labels of the link paths, the segments
#'   that some end paths share; a label holds the labels of the paths it
#'   links, separated by `link_sep`.
#' @param link_sep The separator in the labels of the link paths, matched as
#'   fixed text when the label contains it and as a regular expression
#'   otherwise, so `"|"` and `"\\|"` both work.
#'
#' @return The object with the bins of every end path, from the root on,
#' and their pseudo-bulk counts in the `Dense` slot. The `Design`,
#' `Profile`, `Estimate` and `Significant` slots, computed from
#' earlier bins, are emptied.
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
    .check_strings(
        assay_name = assay_name, root_node = root_node, link_sep = link_sep
    )
    .check(nzchar(link_sep), msg = "'link_sep' must not be an empty string.")
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

    # The labels of a link path, split on link_sep as fixed text when the
    # label contains it, so that "|" works, else as a regular expression,
    # as "\\|" needs
    split_link <- function(link) {
        fixed_sep <- grepl(link_sep, link, fixed = TRUE)
        unlist(strsplit(link, link_sep, fixed = fixed_sep))
    }

    # Every node is a path of the bins, and a link path names two or more
    # of them
    paths <- unique(bin_info[[scmpObj@Parameters@path_col]])
    ends <- unlist(end_node_list)
    links <- unlist(link_node_list)
    .check(length(ends) >= 1L,
        msg = "'end_node_list' must name at least one end path."
    )
    unknown <- setdiff(c(root_node, ends, links), paths)
    .check(length(unknown) == 0L,
        msg = paste0(
            "Not paths of the bins: ", paste(unknown, collapse = ", "),
            ". The paths are: ", paste(paths, collapse = ", "), "."
        )
    )
    for (link in links) {
        parts <- split_link(link)
        .check(length(parts) >= 2L && all(parts %in% paths),
            msg = paste0(
                "The link path '", link, "' must hold the labels of the ",
                "paths it links, separated by link_sep = '", link_sep, "'."
            )
        )
    }
    # Each path is named once: the root among the end paths, an end path
    # named twice or a link path also named as an end path got their bins
    # twice
    named <- c(root_node, ends, links)
    twice <- unique(named[duplicated(named)])
    .check(length(twice) == 0L,
        msg = paste0(
            "Each path must be named once in root_node, end_node_list and ",
            "link_node_list; named more than once: ",
            paste(twice, collapse = ", "), "."
        )
    )
    # The bins of a path named nowhere, such as a linked path left out of
    # end_node_list, are left out
    left_out <- setdiff(paths, named)
    if (verbose && length(left_out) > 0L) {
        message(
            "The bins of paths not named in root_node, end_node_list or ",
            "link_node_list are left out: ", paste(left_out, collapse = ", ")
        )
    }

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

        link_vec <- split_link(link_value)

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
    scmpObj <- .clear_after(scmpObj, "bins")

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
