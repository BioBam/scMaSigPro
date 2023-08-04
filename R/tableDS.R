#' Identification of Major and minor Isoforms in the clusters
#'
#' \code{tableDS} identifies for each Differentialy Spliced Gene (DSG) the clusters where their isoforms belong to, labelling gene transcripts as major (or most expressed) and minor.
#'
#' @title Identification of Major and minor Isoforms in the clusters
#'
#' @description
#' \code{tableDS} identifies for each Differentialy Spliced Gene (DSG) the clusters where their isoforms belong to, labelling gene transcripts as major (or most expressed) and minor.
#'
#' @param seeDS a \code{seeDS} object
#'
#' @details
#' This table includes DSG with 2 or more Isoforms. Mono isoform genes are useful to determine the trends of the cluster. However, as they have only one Isoform, there is not the possibility of comparing minor and major DETs.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{IsoTable}{A classification table that indicates the distribution of isoforms across different clusters.}
#'     \item{IsoClusters}{A data.frame with genes in rows and two columns: first indicates the number of cluster of the major isoform and second the number(s) of cluster(s) of the minor isoforms.}
#'   }
#'
#' @references
#' \itemize{
#'   \item{Nueda, M.J., Martorell, J., Marti, C., Tarazona, S., Conesa, A. 2018. Identification and visualization of differential isoform expression in RNA-seq time series. Bioinformatics. 34, 3, 524-526.}
#'   \item{Nueda, M.J., Tarazona, S., Conesa, A. 2014. Next maSigPro: updating maSigPro bioconductor package for RNA-seq time series. Bioinformatics, 30, 2598-602.}
#'   \item{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006. maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments. Bioinformatics 22, 1096-1102.}
#' }
#'
#' @author Maria Jose Nueda (mj.nueda@ua.es)
#'
#' @seealso \code{\link{seeDS}}, \code{\link{IsoModel}}
#'
#' @examples
#' data(ISOdata)
#' data(ISOdesign)
#' mdis <- make.design.matrix(ISOdesign)
#' MyIso <- IsoModel(data = ISOdata[, -1], gen = ISOdata[, 1], design = mdis, counts = TRUE)
#' Myget <- getDS(MyIso)
#' see <- seeDS(Myget, cluster.all = FALSE, k = 6)
#' table <- tableDS(see)
#' table$IsoTable
tableDS <- function(seeDS) {
  # Extract necessary data from the seeDS object
  Model <- seeDS$Model
  get2 <- seeDS$get2
  cut <- seeDS$cut

  # Extract data, gen, and design from Model
  data <- Model$data
  gen <- Model$gen
  design <- Model$design

  # Extract summary and NumIso.by.gene from get2
  sig.iso2 <- get2$summary
  gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
  NT2 <- seeDS$NumIso.by.gene

  #------------------------------------------------------------
  # Prepare data for clustering
  #------------------------------------------------------------

  # Extract data for significant genes and remove mono-DSG
  data.clust <- get2$sig.genes$sig.profiles
  genes.2 <- names(NT2[NT2 > 1])
  data.clust <- data.clust[gen.sig.iso2 %in% genes.2, ]
  gen.sig.iso22 <- gen.sig.iso2[gen.sig.iso2 %in% genes.2]
  cut <- cut[sig.iso2[gen.sig.iso2 %in% genes.2]]

  # renaming gen.sig.iso22
  gen.sig.iso2 <- gen.sig.iso22

  #-------------------------------------------------------
  # Identifying Major Isoform
  #-------------------------------------------------------

  # Extract unique gene names and their count (G2)
  unic <- unique(gen.sig.iso2)
  G2 <- length(unic)
  Mayor <- NULL
  for (j in 1:G2) {
    # Subset data for the current gene
    zz <- data.clust[gen.sig.iso2 == unic[j], ]
    # Determine the major isoform for the gene and assign names
    Mayor.j <- MayorIso(zz)
    names(Mayor.j) <- rownames(zz)
    Mayor <- c(Mayor, Mayor.j)
  }

  # Prepare data for IsoTable
  cuts <- NULL
  for (i in 1:length(unic)) {
    # Extract the cluster number of the major isoform (cutMi)
    cutMi <- cut[gen.sig.iso2 == unic[i] & Mayor == 1]
    # Extract unique cluster numbers of the minor isoforms (cutmi)
    cutmi <- sort(unique(cut[gen.sig.iso2 == unic[i] & Mayor == 0]))
    if (length(cutmi) > 1) {
      # If there are multiple clusters for minor isoforms, merge them using f3 function
      cutmi <- f3(cutmi)
    }
    # Combine cluster numbers for major and minor isoforms for the gene (cuti)
    cuti <- c(cutMi, cutmi)
    # Append the result to cuts and assign row names
    cuts <- rbind(cuts, cuti)
    rownames(cuts)[i] <- unic[i]
  }

  # Create IsoTable as a data frame
  cuts <- as.data.frame(cuts)
  colnames(cuts) <- c("Cluster.Mayor", "Cluster.minor")
  IsoTable <- table(cuts)

  ## ----------------- RESULTS--------------------------------------------

  # Prepare the output as a list
  out <- list(IsoTable, cuts)
  names(out) <- c("IsoTable", "IsoClusters")
  out
}

#---------------------------------------------------------------------------------------------------
# Auxiliar internal functions: f3, MayorIso
#---------------------------------------------------------------------------------------------------
f3 <- function(x) {
  x <- x[!is.na(x)]
  y <- x[1]
  for (i in 2:length(x)) {
    y <- paste(y, x[i], sep = "&")
  }
  y
}
# ejemplo:
# x<-c(1,2,4,NA,NA)
# f3(x)

#-------------------------------------------------------
MayorIso <- function(zz) {
  if (is.null(nrow(zz))) {
    sol <- 1
  } else {
    M <- apply(zz, 1, sum)
    sol <- as.numeric(M == max(M))
  }
  sol
}
