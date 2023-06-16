make_pseudobulk_design <- function(design.file, paths.vector,
                                   binnedCol = "binnedTime") {
  # Initialize empty df
  pseudo.bulk.profile <- data.frame(NULL)

  # Create Profile
  for (i in paths.vector) {
    # Create a column
    design.file$cell <- rownames(design.file)

    # Create a temp.data.frame per path
    tmp.df <- design.file[design.file[[i]] == 1, ]

    # Select the columns of interest and order
    tmp.df <- tmp.df[order(tmp.df[, binnedCol]), c(binnedCol, "cell")]

    # Group-by
    tmp.df <- tmp.df %>%
      group_by_at(binnedCol) %>%
      summarise(cluster.members = paste0(cell, collapse = "|"))

    tmp.df$bin <- paste0(i, "_bin_", seq(1, nrow(tmp.df)))

    # Add Path Info
    tmp.df$path <- i

    # Add to the frame
    pseudo.bulk.profile <- rbind(pseudo.bulk.profile, tmp.df)
  }

  # Add the cluster size with apply
  pseudo.bulk.profile$bin.size <- apply(pseudo.bulk.profile, 1, calc_bin_size)

  # Add Dummies
  for (i in paths.vector) {
    pseudo.bulk.profile[[i]] <- ifelse(pseudo.bulk.profile$path %in% i, 1, 0)
  }

  # COnvert to a dataframe
  pseudo.bulk.profile <- as.data.frame(pseudo.bulk.profile)

  # Add rownames
  rownames(pseudo.bulk.profile) <- pseudo.bulk.profile$bin

  # Pathway infor
  return(pseudo.bulk.profile)
}
