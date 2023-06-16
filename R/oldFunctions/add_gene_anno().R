add_gene_anno <- function(sim.sce, num.path = 2,
                          base.str = "baseGeneMean",
                          defac.str = "DEFac") {
  # Extract all the gene metadata as a dataframe
  gene.info <- as.data.frame(rowData(sim.sce))

  if (num.path == 2) {
    # Get the columns having information about gene regulation status
    sel.cols <- grep(pattern = paste0(base.str, "|", defac.str), ignore.case = T, value = T, x = colnames(gene.info))

    # Subset
    gene.info <- gene.info[, colnames(gene.info) %in% sel.cols]

    # Order the dataset
    # gene.info <- gene.info[order(
    #     gene.info[[grep(pattern = base.str, ignore.case = T, value = T, x = colnames(gene.info))]], decreasing = T),]
    #
    # Add Additional Columns
    gene.info$PathAExp <- gene.info$BaseGeneMean * gene.info$DEFacPath1
    gene.info$PathBExp <- gene.info$BaseGeneMean * gene.info$DEFacPath2

    # Calculate Change
    gene.info$BaseToPathAChange <- gene.info$PathAExp - gene.info$BaseGeneMean
    gene.info$BaseToPathBChange <- gene.info$PathBExp - gene.info$BaseGeneMean

    # Add Symbols for No-Change
    gene.info[(gene.info$BaseToPathAChange == 0 & gene.info$BaseToPathBChange == 0), "status"] <- "No_Change"

    # Add Changes in one
    gene.info[(gene.info$BaseToPathAChange != 0 & gene.info$BaseToPathBChange == 0), "status"] <- "One_Change"
    gene.info[(gene.info$BaseToPathAChange == 0 & gene.info$BaseToPathBChange != 0), "status"] <- "One_Change"

    # Add Similar Change
    gene.info[(gene.info$BaseToPathAChange < 0 & gene.info$BaseToPathBChange < 0), "status"] <- "Similar_Change"
    gene.info[(gene.info$BaseToPathAChange > 0 & gene.info$BaseToPathBChange > 0), "status"] <- "Similar_Change"

    # Add Opposite Change
    gene.info[(gene.info$BaseToPathAChange < 0 & gene.info$BaseToPathBChange > 0), "status"] <- "Opposite_Change"
    gene.info[(gene.info$BaseToPathAChange > 0 & gene.info$BaseToPathBChange < 0), "status"] <- "Opposite_Change"

    # Add extra column for monocle3
    gene.info$gene_short_name <- rownames(gene.info)

    # Remove if foldChange is less than 1
    gene.info.no.change <- gene.info[gene.info$status == "No_Change", ]
    gene.info.change <- gene.info[!(rownames(gene.info) %in% rownames(gene.info.no.change)), ]
    gene.info.hFC <- gene.info.change[(abs(gene.info.change$BaseToPathAChange) >= 1 | abs(gene.info.change$BaseToPathBChange) >= 1), ]
    gene.info.lFC <- gene.info.change[!(rownames(gene.info.change) %in% rownames(gene.info.hFC)), ]

    # Add another column for status
    gene.info.no.change$status2 <- "No_Change"
    gene.info.hFC$status2 <- "High_FC"
    gene.info.lFC$status2 <- "Low_FC"

    gene.info <- rbind(gene.info.no.change, gene.info.hFC, gene.info.lFC)

    # Add back to the sce object
    return(gene.info)
  }
}
