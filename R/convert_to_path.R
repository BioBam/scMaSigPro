#' Convert Vector Elements to Path Names (Internal)
#'
#' This function transforms a vector by renaming its unique elements
#' (excluding the "root" element) to a sequence named "Path1", "Path2", etc.
#'
#' @param vec A character vector where elements may be repeated and
#' might contain the value "root".
#' @param path_prefix Prefix used to annoate the paths, default is "Path".
#'
#' @return A character vector with the same length as the input where
#' unique elements, excluding "root", are renamed to "Path1", "Path2", etc.
#'
#'
#' @keywords internal
convert_to_path <- function(vec, path_prefix, root_label) {
  # Exclude "root" from the transformation and get the unique values
  unique_vals <- unique(vec[vec != root_label])

  # Create a named vector for the mapping
  name_map <- setNames(paste0(path_prefix, 1:length(unique_vals)), unique_vals)

  # Map and replace the non-root elements
  vec[vec != root_label] <- name_map[vec[vec != root_label]]

  return(vec)
}
