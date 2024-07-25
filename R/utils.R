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


#' Convert List of Named Vectors to Data Frame with Row Names for UpSetR
#'
#' This function converts a list of named vectors to a data frame compatible with UpSetR,
#' retaining the row names corresponding to the unique elements (genes) present in the input list.
#'
#' @param input A list of named vectors to be converted to a data frame compatible with UpSetR.
#'
#' @return A data frame with binary values indicating the presence (1) or absence (0) of each element
#' in the sets, and row names corresponding to the unique elements.
#'
#' @keywords internal
fromListWithNames <- function(input) {
  # Get unique elements (genes)
  elements <- unique(unlist(input))

  # Create binary matrix
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = FALSE))

  # Filter out rows with no data
  data <- data[which(rowSums(data) != 0), ]

  # Set column names
  names(data) <- names(input)

  # Set row names
  rownames(data) <- elements[which(rowSums(data) != 0)]

  return(data)
}
