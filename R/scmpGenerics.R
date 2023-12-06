#' @export
cSparse <- function(object, value = "missing") {
    standardGeneric("cSparse")
}

#' @export
`cSparse<-` <- function(object, value) {
    standardGeneric("cSparse<-")
}

##############################################################################
setMethod(
    "show",
    "scmp",
    function(object) {
        .scmp_show(object)
    }
)

setMethod("cSparse", "scmp", function(object, value) {
    if (identical(value, "missing")) {
        return(object@sparse@colData)  # Getter
    } else {
        object@sparse@colData <- value  # Setter
        return(invisible(object))
    }
})

setReplaceMethod("cSparse", "scmp", function(object, value) {
    object@sparse@colData <- value
    return(object)
})

###############################################################################
scmp <- function(sparse = new("SingleCellExperiment"), 
                 profile = new("sigProfileClass"),
                 estimate = new("estimateClass"),
                 dense = new("SingleCellExperiment"),
                 design = new("designClass"),
                 param = new("paramClass"),
                 sig.genes = new("sigClass")) {
    new("scmp",
        sparse = sparse,
        profile = profile,
        estimate = estimate,
        dense = dense,
        design = design,
        paramClass = param,
        sig.genes = sig.genes # new("sigClass"),
    )
}

