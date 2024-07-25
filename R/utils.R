#' @title Get Identifying the OS from R.
#'
#' @description
#' This function returns the name of the operating system the R session is
#' running on. It differentiates between 'osx', 'linux', and other systems.
#' Visit https://www.r-bloggers.com/
#'
#' @author Will
#'
#' @return A string indicating the operating system: either 'osx', 'linux', or a
#' lower-case version of `Sys.info()['sysname']` if the system is neither OS X
#' nor Linux.
#'
#' @keywords internal
get_os <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") {
      os <- "osx"
    }
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)) {
      os <- "osx"
    }
    if (grepl("linux-gnu", R.version$os)) {
      os <- "linux"
    }
  }
  as.vector(tolower(os))
}

#' Clean or Check Backticks in a String
#'
#' This function either checks for the presence of backticks in a string or removes them, based on the specified action.
#'
#' @param input_string A character string to be processed.
#' @param action A character string specifying the action to perform. Either "check" to check for backticks or "remove" to remove backticks. Defaults to "check".
#'
#' @return If action is "check", returns a logical value indicating the presence of backticks. If action is "remove", returns the input string with all backticks removed.
#'
#' @keywords internal
clean_string <- function(input_string, action = c("check", "remove")) {
  action <- match.arg(action)

  if (action == "check") {
    return(stringr::str_detect(input_string, "`"))
  } else if (action == "remove") {
    # Remove Backticks
    processed_string <- stringr::str_replace_all(input_string, "`", "")
    return(processed_string)
  }
}
