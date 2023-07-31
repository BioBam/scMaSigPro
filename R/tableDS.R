tableDS <- function(seeDS) {
  Model <- seeDS$Model
  get2 <- seeDS$get2
  cut <- seeDS$cut

  data <- Model$data
  gen <- Model$gen
  design <- Model$design

  sig.iso2 <- get2$summary
  gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
  NT2 <- seeDS$NumIso.by.gene

  #------------------------------------------------------------
  data.clust <- get2$sig.genes$sig.profiles

  # Removing mono-DSG:

  genes.2 <- names(NT2[NT2 > 1])
  data.clust <- data.clust[gen.sig.iso2 %in% genes.2, ]
  gen.sig.iso22 <- gen.sig.iso2[gen.sig.iso2 %in% genes.2]
  cut <- cut[sig.iso2[gen.sig.iso2 %in% genes.2]]
  # renaming gen.sig.iso22
  gen.sig.iso2 <- gen.sig.iso22

  #-------------------------------------------------------
  # Identifying Major Isoform
  #-------------------------------------------------------
  unic <- unique(gen.sig.iso2)
  G2 <- length(unic)
  Mayor <- NULL
  for (j in 1:G2)
  {
    zz <- data.clust[gen.sig.iso2 == unic[j], ]
    Mayor.j <- MayorIso(zz)
    names(Mayor.j) <- rownames(zz)
    Mayor <- c(Mayor, Mayor.j)
  }

  cuts <- NULL
  for (i in 1:length(unic))
  {
    cutMi <- cut[gen.sig.iso2 == unic[i] & Mayor == 1]
    cutmi <- sort(unique(cut[gen.sig.iso2 == unic[i] & Mayor == 0]))
    if (length(cutmi) > 1) {
      cutmi <- f3(cutmi)
    }
    cuti <- c(cutMi, cutmi)
    cuts <- rbind(cuts, cuti)
    rownames(cuts)[i] <- unic[i]
  }
  cuts <- as.data.frame(cuts)
  colnames(cuts) <- c("Cluster.Mayor", "Cluster.minor")
  IsoTable <- table(cuts)

  ## ----------------- RESULTS--------------------------------------------

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
  for (i in 2:length(x))
  {
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
