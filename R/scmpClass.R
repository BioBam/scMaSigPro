#' #' scMaSigProClass
#' #'
#' #' A class to represent the ScMaSigPro analysis results and associated data.
#' #' Inherits from \code{SingleCellExperiment}.
#' #'
#' #' @name scMaSigProClass
#' #' @aliases scMaSigProClass
#' #' @slot covariate A \code{covariateClass} object containing covariate information.
#' #' @slot parameters A \code{parameterClass} object containing analysis parameters.
#' #' @slot pVector A \code{pVectorClass} object containing p-vector information.
#' #' @slot tFit A \code{tFitClass} object containing model fitting information.
#' #' @slot sigGenes A \code{sigGeneClass} object containing significant genes information.
#' #' @importFrom S4Vectors DataFrame
#' #' @importFrom SummarizedExperiment SummarizedExperiment
#' #' @importFrom methods setClass
#' #' @importFrom methods setMethod
#' #' @exportClass scMaSigProClass
#'
#' # scMaSigProClass
#' setClass("scMaSigProClass",
#'   contains = c("SingleCellExperiment"),
#'   slots = c(
#'     covariate = "covariateClass",
#'     parameters = "parameterClass",
#'     pVector = "pVectorClass",
#'     tFit = "tFitClass",
#'     sigGenes = "sigGeneClass"
#'   )
#' )
