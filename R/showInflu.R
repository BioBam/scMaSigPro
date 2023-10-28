#' Show or Return the Solution
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @examples
#' \dontrun{
#' # Assuming 'scmpObj' is an object of class 'scMaSigProClass'
#' # with a computed solution:
#' showSol(scmpObj, view = TRUE, return = FALSE)
#' }
#' @export
showInflu <- function(scmpObj, view = TRUE, return = FALSE) {
    # Check Object Validity
    assert_that(is(scmpObj, "scMaSigProClass"),
                msg = "Please provide object of class 'scMaSigPro'"
    )
    
    # Check if the sol exist
    assert_that(!all(dim(scmpObj@scTFit@influ.info) == c(0, 0)),
                msg = "tscore is not computed yet"
    )
    
    # Extract
    influ <- scmpObj@scTFit@influ.info
    
    # If viewing is requested
    if (view) {
        View(influ)
    }
    
    # If requested
    if (return) {
        return(influ)
    }
}
