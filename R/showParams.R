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
#' @importFrom methods slot slotNames
#' @export
showParams <- function(scmpObj, view = FALSE, return = TRUE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Get slot names, assuming 'addParams' is a slot within 'scmpObj'
  all_slots <- slotNames(scmpObj)

  # Get 'addParams' slot data using the correct S4 accessor method
  addParamsData <- slot(scmpObj, "addParams")

  # Get all slots of 'addParams', assuming 'addParams' itself is an S4 object with slots
  params <- lapply(slotNames(addParamsData), function(parameter) {
    slot(addParamsData, parameter)
  })

  # Get the data
  params <- data.frame(
    parameters = slotNames(addParamsData),
    value = unlist(params)
  )

  # If viewing is requested and the 'View' function is available
  if (view && exists("View")) {
    View(params)
  }

  # If requested, return the parameters
  if (return) {
    return(params)
  }
}

# You might need to load necessary packages or define any additional required functions
# like 'assert_that' if not already defined.
