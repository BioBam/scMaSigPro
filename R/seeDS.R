seeDS <- function(get, rsq = 0.4, cluster.all = TRUE, plot.mDSG = FALSE, k = 6,
                  cluster.method = "hclust", k.mclust = FALSE, ...) {
  # cluster.all=TRUE, cluster with significant isoforms of all the genome.
  # cluster.all=FALSE, cluster with DETs of DSGs

  Model <- get$Model
  get2 <- get$get2

  data <- Model$data
  gen <- Model$gen
  design <- Model$design

  sig.iso2 <- get2$summary
  gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
  NT2 <- get$NumIso.by.gene

  #---------------------------------------------------------------------------------------------------
  # cluster.all = TRUE. First p.vector to all the transcripts
  #---------------------------------------------------------------------------------------------------
  if (cluster.all) {
    step3 <- p.vector(data, design, family = Model$pvector2$family)
    Tfit3 <- T.fit(step3)
    get3 <- get.siggenes(Tfit3, vars = "all", rsq = rsq)
    sig.iso3 <- get3$summary

    H <- see.genes(get3$sig.genes, item = "Isoforms", k = k, ...)
    cut <- H$cut[sig.iso2] # tomamos s?lo el cut de las isoformas que nos interesan: 325
  }

  if (!cluster.all) {
    H <- see.genes(get2$sig.genes, item = "Isoforms", k = k, ...)
    cut <- H$cut
  }

  #---------------------------------------------------------------------------------------------------
  # If a plot of mDSG is asked:
  #---------------------------------------------------------------------------------------------------
  if (plot.mDSG) {
    data.clust <- get2$sig.genes$sig.profiles
    genes.1 <- names(NT2[NT2 == 1])
    data.clust1 <- data.clust[gen.sig.iso2 %in% genes.1, ]
    H1 <- see.genes(data.clust1, edesign = design$edesign, cluster.method = cluster.method, k.mclust = k.mclust, k = k, item = "Isoforms", ...)
  }

  ## ----------------- RESULTS--------------------------------------------

  out <- list(Model, get2, NT2, cut, gen.sig.iso2)
  names(out) <- c("Model", "get2", "NumIso.by.gene", "cut", "names.genes")
  out
}
