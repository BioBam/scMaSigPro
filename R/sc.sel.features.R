#' pvector_scMASigPro
#'
#' @name scMaSigPro_t_fit
#' @aliases scMaSigPro_t_fit
#'
#' @param scMaSigPro.obj Description of Slot
#' @param Q Description of Slot
#' @param dist.family Description of Slot
#' @param assay.name Description of Slot
#' @param ep Epsilon
#' @param method Step-bacl
#' @param covar.Q description
#' @param nb.k Dispersion factor negative binomial
#'
# Define a function named scMaSigPro_t_fit
# The function takes in several parameters which includes a single-cell MaSigPro object, a Q value (default is 0.05),
# the distribution family (default is "nb"), assay name (default is "counts"), an epsilon value (default is 0.00001),
# the method (default is "stepback"), and a negative binomial dispersion factor (default is 1)

sc.sel.features <- function(scMaSigPro.obj,
                            Q = 0.05,
                            filter = "rsq",
                            dist.family = "nb",
                            assay.name = "counts",
                            ep = 0.00001,
                            method = "stepback",
                            nb.k = 1) {
  # Vector of accepted distribution families
  filter.accpt <- c("rsq", "pvalue")

  # Set seed for reproducibility
  set.seed(2022)

  # Get the dataframe of of solution
  sol.frame <- showSol(scMaSigPro.obj, view = F)

  # Select the significant models
  if (filter == "rsq") {
    sol.frame.sig <- sol.frame[sol.frame[, 2] >= Q, , drop = F]
  } else if (filter == "pvalue") {
    sol.frame.sig <- sol.frame[sol.frame[, 1] <= Q, , drop = F]
  }

  # Remove P-Value and Models
  sol.frame.sig <- sol.frame.sig[, -c(1:3), drop = F]

  # Add Group Vector Names
  sol.frame.sig.group <- sol.frame.sig

  # Rest the column names
  colnames(sol.frame.sig.group) <- scMaSigPro.obj@covariate@group.vector

  # Get the path specific gene significance
  sel.genes <- sapply(unique(colnames(sol.frame.sig.group)), simplify = T, function(grp.i, grp.pvalue.frame = sol.frame.sig.group) {
    # Get the frame specific p-values
    grp.frame <- grp.pvalue.frame[, colnames(grp.pvalue.frame) == grp.i, drop = F]

    # Select the genes
    sel.genes <- rownames(grp.frame[which(apply(grp.frame, 1, function(x) {
      any(!is.na(x))
    })), ])

    # Return the gene vector
    return(sel.genes)
  })

  sel.genes.frame <- as.data.frame(sel.genes)

  # Path2 and Path3vsPath2
  sig.gene.list <- list()

  sig.gene.list[[colnames(sel.genes)[1]]] <- setdiff(sel.genes[, 1], sel.genes[, 2])
  sig.gene.list[[colnames(sel.genes)[2]]] <- setdiff(sel.genes[, 2], sel.genes[, 1])
  sig.gene.list[["shared"]] <- intersect(sel.genes[, 2], sel.genes[, 1])

  # Create Tfit Object
  sigGene.obj <- new("sigGeneClass",
    group.associations = sig.gene.list,
    sel.genes.frame = sel.genes.frame
  )

  # Add to MaSigPro Object
  scMaSigPro.obj@sigGenes <- sigGene.obj

  return(scMaSigPro.obj)
}
