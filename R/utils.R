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
