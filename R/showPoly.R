#' Show the terms of the polynomial term
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#'
#' @return Return the terms of the polynomial model.
#'
#' @export
showPoly <- function(scmpObj) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(all(!is.na(colnames(scmpObj@edesign@dis)) | length(colnames(scmpObj@edesign@dis) > 1)),
    msg = "Please setup the model first, using 'sc.make.design.matrix()'"
  )

  # Extract columns
  df.col <- colnames(scmpObj@edesign@dis)

  # Extract betas
  beta_names <- paste0("beta", seq(1:length(df.col)))

  # Generate formula string
  formula_parts <- vapply(seq_along(df.col),
    function(i) paste(beta_names[i], "*", df.col[i], sep = ""),
    FUN.VALUE = character(1)
  )

  # Make formula
  formula_string <- paste("beta0", paste(formula_parts, collapse = " + "), sep = " + ")

  return(formula_string)
}
