#' Select the Longer of Two Vectors
#'
#' This function compares the lengths of two vectors and returns the longer one.
#' If both vectors have the same length, the function returns 1.
#'
#' @param vector1 A numeric vector.
#' @param vector2 A numeric vector.
#' @param vector1_label A label (character string) for `vector1`.
#' @param vector2_label A label (character string) for `vector2`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{long_vec}: The longer vector of `vector1` and `vector2`.
#'   \item \code{long_vec_label}: The label corresponding to the longer vector.
#'   \item \code{short_vec}: The shorter vector of `vector1` and `vector2`.
#'   \item \code{short_vec_label}: The label corresponding to the shorter vector.
#' }
#' If both vectors are of the same length, the function returns 1.
#'
#'


select_longer_vector <- function(vector1, vector2,
                                 vector1_label, vector2_label) {
    if (length(vector1) > length(vector2)) {
        return(list(long_vec = vector1, long_vec_label = vector1_label,
                    short_vec = vector2, short_vec_label = vector2_label
        ))
    } else if (length(vector2) > length(vector1)) {
        return(list(long_vec = vector2, long_vec_label = vector2_label,
                    short_vec = vector1, short_vec_label = vector1_label
        ))
    } else {
        return(1)
    }
}
