#' scMaSigProClass
#' @name scMaSigProClass
#' @aliases scMaSigProClass
#' @slot covariate Description of Slot.
#' @slot parameters Description of Slot.
#' @exportClass scMaSigProClass
#'

# Object for MaSigPro
setClass("covariateClass",
  slots = c(
    factor.design = "matrix",
    time.series = "numeric",
    covariate = "matrix",
    covariate.vector = "character",
    group.vector = "character",
    groups = "character"
  )
)

# scMaSigProParameterClass
setClass("parameterClass",
  slots = c(
    path.col = "character",
    time.col = "character",
    poly.order = "integer",
    assay.name = "character",
    p.vector.sig = "numeric",
    m.t.c.method = "character",
    min.counts = "numeric",
    model.type = "character",
    model.control = "list",
    covar.sig = "numeric"
  )
)

# Define the constructor method
setMethod("initialize", "parameterClass", function(.Object, path.col, time.col = "Pseudotime",
                                                   poly.order = 2, assay.name = "counts",
                                                   p.vector.sig = 0.05,
                                                   m.t.c.method = "BH",
                                                   min.counts = 6, model.type,
                                                   model.control,
                                                   covar.sig = 0.05) {
  .Object@path.col <- path.col
  .Object@time.col <- time.col
  .Object@poly.order <- poly.order
  .Object@assay.name <- assay.name
  .Object@p.vector.sig <- p.vector.sig
  .Object@m.t.c.method <- m.t.c.method
  .Object@min.counts <- min.counts
  .Object@covar.sig <- covar.sig
  return(.Object)
})


# pVectorClass
setClass("pVectorClass",
  slots = c(
    full.models = "list",
    intercept.models = "list",
    p.value = "numeric",
    adj.p.value = "numeric"
  )
)

# sigGeneClass
setClass("sigGeneClass",
  slots = c(
    group.associations = "list",
    sel.genes.frame = "list"
  )
)

# T.fitClass
setClass("tFitClass",
  slots = c(
    covar.models = "list",
    model.attributes = "list"
  )
)

# scMaSigProClass
setClass("scMaSigProClass",
  contains = c("SingleCellExperiment"),
  slots = c(
    covariate = "covariateClass",
    parameters = "parameterClass",
    pVector = "pVectorClass",
    tFit = "tFitClass",
    sigGenes = "sigGeneClass"
  )
)
