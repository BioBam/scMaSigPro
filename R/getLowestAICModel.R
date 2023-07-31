#' Get Model with Lowest AIC Function
#'
#' Description:
#' This function takes a list of fitted models and returns the name of the model
#' with the lowest AIC (Akaike Information Criterion) value.
#'
#' Usage:
#' getLowestAICModel(model_list)
#'
#' Arguments:
#'   \describe{
#'     \item{model_list}{A list of fitted models for which you want to find the
#'                      one with the lowest AIC. Each element of the list should be
#'                      a fitted model object that supports the AIC function.}
#'   }
#'
#' Details:
#' The function computes the AIC values for all the models in the input list
#' using the `AIC` function. It then finds the index of the model with the
#' lowest AIC value and retrieves its name.
#'
#' Value:
#' The function returns the name of the model with the lowest AIC value as a character.
#'
#' @examples
#' # Fit models and create a list
#' model1 <- lm(mpg ~ wt, data = mtcars)
#' model2 <- lm(mpg ~ wt + hp, data = mtcars)
#' model3 <- lm(mpg ~ wt + hp + qsec, data = mtcars)
#' model_list <- list(model1, model2, model3)
#'
#' # Get the model with the lowest AIC
#' lowest_aic_model <- getLowestAICModel(model_list)
#' print(lowest_aic_model)
#'
#' @importFrom stats AIC
#'
#' @keywords internal
getLowestAICModel <- function(model_list) {
  # Compute AIC values for all models
  aic_values <- sapply(model_list, AIC)

  # Find index of model with the lowest AIC
  lowest_aic_index <- which.min(aic_values)

  # Get the name
  lowest_aic_name <- names(model_list)[lowest_aic_index]

  # Return the model with the lowest AIC
  return(lowest_aic_name)
}
