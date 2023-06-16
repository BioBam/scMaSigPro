classify_sig_genes <- function(masigpro.tstep, rsq.value, draw.venn = F,
                               count.table, num_group = 2) {
  # Identify significant genes
  invisible(capture.output(sigs <- get.siggenes(masigpro.tstep,
    rsq = rsq.value,
    vars = "groups", r = 0.7
  )))


  if (is.null(sigs$sig.genes)) {
    return(list(
      detected.genes = "No Significant Genes found",
      undetected.gene = rownames(count.table)
    ))
  } else {
    # Path2 and Path3vsPath2
    shared <- intersect(sigs$summary[, 2], sigs$summary[, 1])

    # Path2
    grp1 <- setdiff(sigs$summary[, 1], sigs$summary[, 2])

    # Path3vsPath2
    grp2 <- setdiff(sigs$summary[, 2], sigs$summary[, 1])

    # All detected
    detected <- c(shared, grp2, grp1)

    # Remove empty element
    detected <- detected[!(detected == " ")]

    # Get the undetected genes
    undetected <- rownames(count.table)[!(rownames(count.table)) %in% detected]

    if (draw.venn == T & length(detected) > 1) {
      # Make Venn Diagram
      maSigPro::suma2Venn(sigs$summary[, c(1:2)])
    }

    # Return
    return(list(detected.genes = detected, undetected.gene = undetected))
  }
}
