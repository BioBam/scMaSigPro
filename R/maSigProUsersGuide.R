maSigProUsersGuide <- function(view = TRUE) {
  f <- system.file("doc", "maSigProUsersGuide.pdf", package = "maSigPro")
  if (view) {
    if (.Platform$OS.type == "windows" && interactive()) {
      shell.exec(f) # consider using browseUrl() instead?
    } else {
      system(paste(Sys.getenv("R_PDFVIEWER"), f, "&"))
    }
  }
  return(f)
}
