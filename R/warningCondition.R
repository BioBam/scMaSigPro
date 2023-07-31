#' Custom Warning Function
#'
#' This function generates a warning with a custom message if a given condition is TRUE.
#'
#' @param condition A logical value indicating whether the warning should be triggered.
#' @param message A character string specifying the warning message.
#'   
#' @examples
#' # Generate a warning if a condition is TRUE
#' warningCondition(TRUE, "This is a custom warning message.")
#' 
#' # Generate no warning if the condition is FALSE
#' warningCondition(FALSE, "This message won't trigger a warning.")
#'
#' @keywords internal
#' @export
warningCondition <- function(condition, message) {
    if (condition) {
        warning(message, call. = FALSE)
    }
}
