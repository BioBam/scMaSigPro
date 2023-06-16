make_bulk_counts <- function(counts, cluster_member_col = "cluster.members",
                             bin_col = "bin", pseudo_bulk_profile, cluster.count.by = "sum") {
  meta.info <- pseudo_bulk_profile[, c(cluster_member_col, bin_col)]

  # initialize an empty dataframe
  pseudo_bulk_counts <- data.frame(matrix(nrow = nrow(counts), ncol = 0))

  # One row at a time
  for (i in c(1:nrow(meta.info))) {
    # Select Rows One by One
    selRow <- meta.info[i, , drop = F]

    # Save the column Names
    col_name <- selRow[, 2, drop = F]

    # cluster components
    clus_com <- c(stringr::str_split(selRow[, 1], "\\|"))[[1]]

    # Select columns from the frame
    sel_cols_df <- counts[, colnames(counts) %in% clus_com, drop = F]

    # Find mean or sum rowise
    if (cluster.count.by == "mean") {
      sel_cols_df_mean <- as.data.frame(rowMeans(sel_cols_df))
    } else if (cluster.count.by == "sum") {
      sel_cols_df_mean <- as.data.frame(rowSums(sel_cols_df))
    } else {
      stop()
    }

    # Put the new column name
    colnames(sel_cols_df_mean) <- col_name

    # Join the frames
    pseudo_bulk_counts <- cbind(pseudo_bulk_counts, sel_cols_df_mean)
  }

  # return
  return(pseudo_bulk_counts)
}
