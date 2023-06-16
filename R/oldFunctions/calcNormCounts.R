# Seurat normCounts
calcNormCounts <- function(rawCounts, cat, size_fac = 10000) {
  # Base
  suppressPackageStartupMessages(require(Seurat))

  # Make Seurat Object
  seuratObject <- CreateSeuratObject(counts = rawCounts, assay = "scRNA")

  ## Lib.size
  if (cat == "libSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject,
      normalization.method = "RC", scale.factor = size_fac
    )
    # Extract Counts
    seuratNormCounts <- seuratNormObject@assays$scRNA@data
    # Return dgCMatrix
    return(seuratNormCounts)
  }

  ## Log + Lib.size
  else if (cat == "logLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "LogNormalize"
    )
    # Extract Counts
    seuratNormCounts <- seuratNormObject@assays$scRNA@data
    # Return dgCMatrix
    return(seuratNormCounts)
  }

  ## Centered + Lib.size
  else if (cat == "cenLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "RC"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$scRNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = F)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "scLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "RC"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$scRNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  }
  ## Centered + log + Lib.size
  else if (cat == "cenLogLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "LogNormalize"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$scRNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = F)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "scLogLibSize") {
    # Normalize
    seuratNormObject <- NormalizeData(
      object = seuratObject, scale.factor = size_fac,
      normalization.method = "LogNormalize"
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$scRNA@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "sctransform") {
    # Call ScTranform
    suppressPackageStartupMessages(require(sctransform))

    # Normalize
    seuratNormObject <- SCTransform(seuratObject,
      assay = "scRNA",
      do.center = F,
      do.scale = F,
      # vars.to.regress = "percent.mt",
      verbose = FALSE
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$SCT@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "cenSctransform") {
    # Call ScTranform
    suppressPackageStartupMessages(require(sctransform))

    # Normalize
    seuratNormObject <- SCTransform(seuratObject,
      assay = "scRNA",
      do.center = T,
      do.scale = F,
      # vars.to.regress = "percent.mt",
      verbose = FALSE
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$SCT@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else if (cat == "scSctransform") {
    # Call ScTranform
    suppressPackageStartupMessages(require(sctransform))

    # Normalize
    seuratNormObject <- SCTransform(seuratObject,
      assay = "scRNA",
      do.center = T,
      do.scale = T,
      # vars.to.regress = "percent.mt",
      verbose = FALSE
    )
    # Extract Counts
    seuratNormCounts <- t(as.matrix(seuratNormObject@assays$SCT@data))
    seuratNormCounts <- scale(seuratNormCounts, center = TRUE, scale = T)
    seuratNormCounts <- t(seuratNormCounts)
    # Return dgCMatrix
    return(seuratNormCounts)
  } else {
    stop()
  }
}
