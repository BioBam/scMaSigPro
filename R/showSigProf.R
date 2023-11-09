#' Show or Return the Solution
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param includeInflu logical, whether to add gene with influential data in the solution.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @examples
#' \donttest{
#' # Assuming 'scmpObj' is an object of class 'scMaSigProClass'
#' # with a computed solution:
#' showSol(scmpObj, view = TRUE, return = FALSE)
#' }
#' @importFrom utils View
#' @export
showSigProf <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = FALSE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@sol) == c(0, 0)),
    msg = "Sol is not computed yet"
  )

  # Extract
  sol <- showSol(scmpObj, view = FALSE, return = TRUE, influ = influ) %>% as.data.frame()
  # Extract rownames
  bulk.counts <- scmpObj@compress.sce@assays@data@listData$bulk.counts
  bulk.counts <- bulk.counts[rownames(bulk.counts) %in% rownames(sol), , drop = FALSE]

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    bulk.counts <- bulk.counts[!(rownames(bulk.counts) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(as.matrix(bulk.counts))
  }

  # If requested
  if (return) {
    return(bulk.counts)
  }
}
