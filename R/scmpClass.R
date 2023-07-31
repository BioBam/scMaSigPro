#' scMaSigProClass
#'
#' A class to represent the ScMaSigPro analysis results and associated data.
#' Inherits from \code{SingleCellExperiment}.
#'
#' @name scMaSigProClass
#' @aliases scMaSigProClass
#' @slot covariate A \code{covariateClass} object containing covariate information.
#' @slot parameters A \code{parameterClass} object containing analysis parameters.
#' @slot pVector A \code{pVectorClass} object containing p-vector information.
#' @slot tFit A \code{tFitClass} object containing model fitting information.
#' @slot sigGenes A \code{sigGeneClass} object containing significant genes information.
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods setClass
#' @importFrom methods setMethod
#' @exportClass scMaSigProClass

# Define covariateClass
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

#' parameterClass
#'
#' A class to represent the analysis parameters for ScMaSigPro.
#'
#' @name parameterClass
#' @slot path.col A character vector specifying the column name for the path data.
#' @slot time.col A character vector specifying the column name for the time data.
#' @slot poly.order An integer specifying the polynomial order for the time covariate.
#' @slot assay.name A character vector specifying the assay name.
#' @slot p.vector.sig A numeric vector specifying the significance level for p-vector.
#' @slot m.t.c.method A character vector specifying the method for multiple testing correction.
#' @slot min.counts A numeric vector specifying the minimum counts for genes.
#' @slot model.type A character vector specifying the model type.
#' @slot model.control A list specifying the model control parameters.
#' @slot covar.sig A numeric vector specifying the significance level for covariate modeling.
#' @importFrom methods setClass
#' @importFrom methods setMethod
#' @exportClass parameterClass

# Define parameterClass
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

# Define the constructor method for parameterClass
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


# Define pVectorClass
setClass("pVectorClass",
  slots = c(
    full.models = "list",
    intercept.models = "list",
    p.value = "numeric",
    adj.p.value = "numeric"
  )
)

#' sigGeneClass
#'
#' A class to represent significant genes and their associations.
#'
#' @name sigGeneClass
#' @slot group.associations A list representing group associations.
#' @slot sel.genes.frame A list representing selected genes frame.
#' @importFrom methods setClass
#' @importFrom methods setMethod
#' @exportClass sigGeneClass

# Define sigGeneClass
setClass("sigGeneClass",
  slots = c(
    group.associations = "list",
    sel.genes.frame = "list"
  )
)

#' tFitClass
#'
#' A class to represent model fitting information.
#'
#' @name tFitClass
#' @slot covar.models A list representing covariate models.
#' @slot model.attributes A list representing model attributes.
#' @importFrom methods setClass
#' @importFrom methods setMethod
#' @exportClass tFitClass

# Define tFitClass
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

#' Create scMaSigProClass Object
#'
#' A constructor function to create an instance of \code{scMaSigProClass}.
#'
#' @param SingleCellExperiment An instance of \code{SingleCellExperiment}.
#' @param covariate An instance of \code{covariateClass}.
#' @param parameters An instance of \code{parameterClass}.
#' @param pVector An instance of \code{pVectorClass}.
#' @param tFit An instance of \code{tFitClass}.
#' @param sigGenes An instance of \code{sigGeneClass}.
#'
#' @return An instance of \code{scMaSigProClass}.
#' @importFrom methods new
#' @export
newScMaSigPro <- function(SingleCellExperiment, covariate, parameters, pVector, tFit, sigGenes) {
  new("scMaSigProClass",
    SingleCellExperiment = SingleCellExperiment,
    covariate = covariate,
    parameters = parameters,
    pVector = pVector,
    tFit = tFit,
    sigGenes = sigGenes
  )
}
