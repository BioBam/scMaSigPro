#' @title Plot Model Diagnostics
#'
#' @description
#' Plots model residuals for the requested gene. Can be used for model
#' diagnostics of optimized model or full model.
#'
#' @importFrom stats residuals fitted rstandard predict
#' @importFrom stringr str_remove
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param feature_id Name of the gene to be plotted.
#' @param model Type of model to be used for plotting. Can be either 'optimized',
#' 'intercept'and 'full'. Default is 'optimized'.
#'
#' @return ggplot2 plot object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
#'
plotDiagnostics <- function(scmpObj,
                            feature_id,
                            model = "optimized") {
  # Check if the results exist
  assert_that(!isEmpty(scmpObj@Significant@genes),
    msg = paste("No Significant genes found, please run the workflow first")
  )

  # Check model type
  assert_that(all(model %in% c("optimized", "full", "intercept")),
    msg = paste("The requested gene is not available in the significant genes list")
  )

  # Get all the available genes
  avail_genes <- unique(unlist(scmpObj@Significant@genes))

  # Check if the requested gene is availble as singificant gene
  assert_that(all(feature_id %in% avail_genes),
    msg = paste("The requested gene is not available in the significant genes list")
  )

  # Get the counts for the genes
  y_mat <- as.matrix(scmpObj@Dense@assays@data$bulk.counts)

  # Extract gene counts
  y_df <- y_mat[feature_id, , drop = TRUE] %>%
    as.data.frame()
  colnames(y_df) <- "y"

  # Set theme for the plot
  theme_set <- theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(
        color = "grey90", linewidth = 0.3, linetype = "dashed"
      ),legend.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  # Based on model
  if (model == "full" || model == "intercept") {
    # Create Model input
    regression_matrix <- scmpObj@Design@predictor_matrix
  } else {
    # Get coefficients
    coeff_matrix <- showCoeff(scmpObj) %>%
      data.frame()

    # Subset
    coeff_matrix_sub <- coeff_matrix[feature_id, , drop = FALSE]
    all_terms <- colnames(coeff_matrix_sub)[-1]

    # Remove any occurance of beta
    all_terms <- str_remove(string = all_terms, pattern = "beta")

    # Remove intercept
    all_terms <- all_terms[-1]

    # Get the terms
    terms <- coeff_matrix_sub[apply(coeff_matrix_sub, 2, FUN = function(i) {
      i != 0
    })]

    # Extract the prediction
    prediction_matrix <- scmpObj@Design@predictor_matrix

    # Create Model input
    regression_matrix <- prediction_matrix[, all_terms, drop = FALSE]
  }

  # Combine regression matrix with counts
  model_df <- cbind(y_df, regression_matrix)

  # Extract offsets
  if (scmpObj@Parameters@offset) {
    offset_data <- scmpObj@Design@offset
  } else {
    offset_data <- NULL
  }

  if (model == "intercept") {
    glm_s3 <- glm(y ~ 1,
      data = as.data.frame(model_df),
      family = scmpObj@Parameters@distribution,
      epsilon = scmpObj@Parameters@epsilon,
      offset = offset_data,
      weights = NULL,
      maxit = scmpObj@Parameters@max_it
    )
  } else {
    # Fit Model
    glm_s3 <- glm(y ~ .,
      data = as.data.frame(model_df),
      family = scmpObj@Parameters@distribution,
      epsilon = scmpObj@Parameters@epsilon,
      offset = offset_data,
      weights = NULL,
      maxit = scmpObj@Parameters@max_it
    )
  }
  # Extract residuals
  raw_residuals <- residuals(glm_s3, type = "response")
  pearson_residuals <- residuals(glm_s3, type = "pearson")
  deviance_residuals <- residuals(glm_s3, type = "deviance")
  standardized_residuals <- rstandard(glm_s3)

  # Extract fitted values
  fitted_values <- fitted(glm_s3)

  # Calculate predicted responses
  predicted_responses <- predict(glm_s3, type = "response")

  # Q-Q plot
  # This plot of the standardized residuals allows you to check if the
  # residuals are approximately normally distributed, which is an
  # Calculate standardized residuals

  # Create a dataframe for plotting
  qq_df <- data.frame(std_residuals = standardized_residuals)

  # Create the Normal Q-Q plot
  qq_plot <- ggplot(qq_df, aes(sample = .data$std_residuals)) +
    stat_qq(color = "#EE446F", size = 2) +
    stat_qq_line(
      color = "#15918A", alpha = 1,
      linewidth = 1, linetype = "solid"
    ) +
    labs(
      subtitle = feature_id,
      title = "A. Normal Q-Q Plot of Standardized Residuals",
      x = "Theoretical Quantiles",
      y = "Standardized Residuals"
    ) +
    theme_set

  # Residuals vs Fitted Plot:
  # This plot helps to detect non-linearity, unequal error variances,
  # and outliers. It is used to check the assumption of homoscedasticity.
  # Assuming your glm model is stored in 'glm_s3'

  # Create a dataframe for plotting
  res_vs_fit_df <- data.frame(
    Fitted = fitted_values,
    StdResiduals = standardized_residuals
  )

  # Create the Residuals vs Fitted plot
  stdRes_vs_fitted_plot <- ggplot(res_vs_fit_df, aes(x = .data$Fitted, y = .data$StdResiduals)) +
    geom_point(color = "#EE446F", size = 2) + # Plot the points
    geom_hline(yintercept = 0, linetype = "dashed", color = "#9F7BB8") + # Add a horizontal line at 0
    labs(
      subtitle = feature_id,
      title = "B. Standardized Residuals vs Fitted Plot", x = "Fitted Values", y = "Standardized Residuals"
    ) +
    theme_set +
    geom_smooth(
      method = "loess", color = "#15918A", alpha = 0.4, formula = y ~ x,
      linewidth = 1, linetype = "solid", se = TRUE, fill = "#FDC659"
    ) # Add a loess smoothed line

  # Scale-Location Plot (or Spread-Location Plot): This shows if residuals
  # are spread equally along the ranges of predictors. This is a way to check
  # the assumption of equal variance (homoscedasticity).
  # Calculate square root of absolute standardized residuals
  sqrt_abs_std_residuals <- sqrt(abs(standardized_residuals))

  # Create a dataframe for plotting
  scale_loc_df <- data.frame(
    Fitted = fitted_values,
    SqrtAbsStdResiduals = sqrt_abs_std_residuals
  )

  # Create the Scale-Location plot
  scale_loc_plot <- ggplot(scale_loc_df, aes(x = .data$Fitted, y = .data$SqrtAbsStdResiduals)) +
    geom_point(color = "#EE446F", size = 2) + # Plot the points
    geom_smooth(
      method = "loess", color = "#15918A", formula = y ~ x,
      linewidth = 1, linetype = "solid", se = TRUE,
      fill = "#FDC659", alpha = 0.4
    ) +
    labs(
      subtitle = feature_id,
      title = "C. Scale-Location Plot",
      x = "Fitted Values",
      y = "Sqrt of Absolute Standardized Residuals"
    ) +
    theme_set

  # A Dfun plot is less standard and typically refers to a plot that can help
  # assess the appropriateness of the link function in a GLM. This might involve
  # plotting observed versus expected responses or derivatives of the link function.
  # The exact implementation can depend on the specific family and link function used
  # in your model. For a basic implementation, you might plot observed vs. predicted responses.

  # Create a dataframe for plotting
  dfun_df <- data.frame(Observed = model_df$y, Predicted = predicted_responses)

  # Create the Dfun plot
  defunc_plot <- ggplot(dfun_df, aes(x = .data$Observed, y = .data$Predicted)) +
    geom_point(color = "#EE446F", size = 2) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "solid", color = "#15918A",
      linewidth = 1, alpha = 1
    ) +
    labs(
      subtitle = feature_id,
      title = "D. Observed vs Predicted",
      x = "Observed",
      y = "Predicted"
    ) +
    theme_set

  # Return Plots
  return((qq_plot + stdRes_vs_fitted_plot) / (scale_loc_plot + defunc_plot))
}
