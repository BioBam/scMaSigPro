#' @title estimateSizeFactorsForMatrix from DESeq2
#' @author Simon Anders
#' Please cite as
#' Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and
#' dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#'
#' Low-level function to estimate size factors with robust regression.
#'
#' @description
#' Given a matrix or data frame of count data, this function estimates the size
#' factors as follows: Each column is divided by the geometric means of the
#' rows. The median (or, if requested, another location estimator) of these
#' ratios (skipping the genes with a geometric mean of zero) is used as the size
#' factor for this column.
#'
#' @param counts a matrix or data frame of counts, i.e., non-negative integer
#' values
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used. However, especially for low counts, the
#' function from `genefilter` package may give better results.
#' @param geoMeans by default this is not provided, and the
#' geometric means of the counts are calculated within the function.
#' A vector of geometric means from another count matrix can be provided
#' for a "frozen" size factor calculation
#' @param controlGenes optional, numeric or logical index vector specifying those genes to
#' use for size factor estimation (e.g. housekeeping or spike-in genes)
#' @param type standard median ratio (\code{"ratio"}) or where the
#' geometric mean is only calculated over positive counts per row
#' (\code{"poscounts"})
#' @return a vector with the estimates size factors, one element per column
#'
#' @importFrom MatrixGenerics rowMeans rowSums
#' @keywords internal
scmp_estimateSizeFactorsForMatrix <- function(counts, locfunc = stats::median,
                                              geoMeans, controlGenes,
                                              type = c("ratio", "poscounts")) {
  type <- match.arg(type, c("ratio", "poscounts"))
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    if (type == "ratio") {
      loggeomeans <- MatrixGenerics::rowMeans(log(counts))
    } else if (type == "poscounts") {
      lc <- log(counts)
      lc[!is.finite(lc)] <- 0
      loggeomeans <- MatrixGenerics::rowMeans(lc)
      allZero <- MatrixGenerics::rowSums(counts) == 0
      loggeomeans[allZero] <- -Inf
    }
  } else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if (!(is.numeric(controlGenes) | is.logical(controlGenes))) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes, , drop = FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  if (incomingGeoMeans) {
    # stabilize size factors to have geometric mean of 1
    sf <- sf / exp(mean(log(sf)))
  }
  sf
}
