PodiumChange <- function(get, only.sig.iso = FALSE, comparison = c("any", "groups", "specific"), group.name = "Ctr", time.points = 0) {
  Model <- get$Model
  get2 <- get$get2

  data <- Model$data
  gen <- Model$gen

  edesign <- Model$design$edesign
  repvect <- edesign[, 2]

  if (only.sig.iso) {
    sig.iso2 <- get2$summary
    gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
    NT2 <- get$NumIso.by.gene
    data.clust <- as.matrix(get2$sig.genes$sig.profiles)
    # Removing monoIsoform Genes (there is not podium change)
    genes.2 <- names(NT2[NT2 > 1])
    data.clust <- data.clust[gen.sig.iso2 %in% genes.2, ]
    gen.sig.iso22 <- gen.sig.iso2[gen.sig.iso2 %in% genes.2]
    # renaming gen.sig.iso22
    gen.sig.iso2 <- gen.sig.iso22
  } else {
    data.clust <- as.matrix(data[gen %in% Model$DSG, ])
    sig.iso2 <- rownames(data.clust)
    gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
    NT2 <- tapply(sig.iso2, gen.sig.iso2, length)
    # Here, there is not any mDSG because in this analysis it is considered only (>1 iso)
  }

  if (length(gen.sig.iso2) == 0) {
    print("No selected genes with more than 1 Isoform")
  } else {
    #-------------------------------------------------------
    # Major Isoform identification
    #-------------------------------------------------------

    time.M <- tapply(edesign[, 1], repvect, mean)
    if (ncol(edesign) == 3) {
      name3 <- colnames(edesign)[3]
      edesign3 <- as.data.frame(edesign[, 3:ncol(edesign)])
      colnames(edesign3) <- name3
    } else {
      edesign3 <- edesign[, 3:ncol(edesign)]
    }
    groups.M <- apply(edesign3, 2, function(x) {
      tapply(x, repvect, mean)
    })

    unic <- unique(gen.sig.iso2)
    Mayor <- NULL
    LIST <- NULL

    for (i in 1:length(unic))
    {
      zz <- data.clust[gen.sig.iso2 == unic[i], ]
      M <- MayorIso(zz)
      zzM <- zz[M == 1, ]
      MzzM <- tapply(zzM, repvect, mean)
      zzm <- zz[M != max(M), ]

      if (is.null(nrow(zzm))) {
        ni <- 1
      } else {
        ni <- nrow(zzm)
      }

      if (ni == 1) {
        Mzzm <- tapply(zzm, repvect, mean)
      } else {
        Mzzm <- t(apply(zzm, 1, function(x) {
          tapply(x, repvect, mean)
        }))
      }

      if (ni == 1) {
        dif <- MzzM - Mzzm
      } else {
        dif <- t(apply(Mzzm, 1, function(x) {
          MzzM - x
        }))
      }
      # Comparison = "any" ----------------------------------------------
      if (comparison == "any") {
        if (any(dif < 0)) {
          LIST <- c(LIST, unic[i])
        }
      }
      # Comparison = "specific" ----------------------------------------------
      else if (comparison == "specific") {
        col <- groups.M[, colnames(groups.M) == group.name]
        if (ni == 1) {
          change <- all(dif[col == 1 & time.M == time.points] < 0)
        } else {
          change <- apply(dif[, col == 1 & time.M == time.points], 1, function(x) {
            all(x < 0)
          })
        }
        if (any(change)) LIST <- c(LIST, unic[i])
      }
      # Comparison = group ----------------------------------------------
      else if (comparison == "group") {
        mayors <- NULL
        for (k in 3:ncol(edesign))
        {
          mayors <- cbind(mayors, MayorIso(zz[, edesign[, k] == 1]))
        }
        # When all columns match, substraction with any of them will be 0:
        if (all(mayors - mayors[, 1] != 0)) LIST <- c(LIST, unic[i])
      }
    }

    # lists of genes and isoforms:
    gen.L <- gen.sig.iso2[gen.sig.iso2 %in% LIST]
    data.L <- data.clust[gen.sig.iso2 %in% LIST, ]

    output <- list(LIST, data.L, gen.L, edesign)
    names(output) <- c("L", "data.L", "gen.L", "edesign")
    output
  }
}
