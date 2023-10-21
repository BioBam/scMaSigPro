#' @title Class "sigClass"
#'
#' @description A class for Storing significant genes
#'
#' @slot summary A dataframe of genes
#' @slot sig.genes A list
#'
#' @name sigClass
#' @aliases sigClass-class
#' @rdname sigClass-class
#' @exportClass sigClass
#' @keywords classes
#'
setClass(
  "sigClass",
  representation(
    summary = "data.frame",
    sig.genes = "list"
  ),
  validity = function(object) {
    if (!is.data.frame(object@summary)) {
      stop("summary slot must be a data.frame")
    }
    if (!is.list(object@sig.genes)) {
      stop("sig.genes slot must be a list")
    }
  },
  prototype = list(
    summary = data.frame(),
    sig.genes = list()
  )
)
