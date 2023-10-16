#' Visually Inspect Model Fit
#'
#' This function is used to visually inspect the model fit of a generalized linear
#' model (GLM) and provides different types of plots for validation purposes.
#'
#' @param glm.model A model from the output of \link[stats]{glm}.
#' @param plot.type Plot type. Currently supported options are:
#'   - "fitted_vs_pearson": Fitted values vs Pearson residuals plot.
#'   - "cooks_vs_obs": Cook's distance vs observations bar plot.
#'   - "pearson_vs_obs": Pearson residuals vs observed values plot.
#'   - "fitted_vs_obs": Fitted values vs observed values plot.
#'   - "std_vs_fitted": Standard residuals vs fitted values plot.
#' @param point.color Point colour for the plots.
#' @param point.shape Point shape for the plots.
#' @param point.alpha Point alpha for the plots.
#' @param ref.line.y.intercept Y intercept of the reference line in some plots.
#' @param func Which function is used (options: "glm" or "zeroinfl").
#' @param input.data No information.
#' @param bar_width Bar width for the Cook's distance plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_bar facet_wrap
#' @importFrom stats fitted resid rstandard
#'
#' @keywords internal
#' @export
glm_inspection <- function(glm.model, plot.type = "fitted_vs_pearson",
                           point.color = "blue", point.shape = 3, bar_width = 0.05,
                           point.alpha = 0.8, ref.line.y.intercept = 0,
                           func = "glm", input.data = NULL) {
  # Create a base Plot
  base_plot <- ggplot() +
    theme_classic() +
    theme(
      panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = rel(0.7)),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )

  if (plot.type == "fitted_vs_pearson") {
    # Extract Fitted Value
    fit_vals <- fitted(glm.model)

    # Extract pearson residuals
    pear_resid <- resid(glm.model, type = "pearson")

    # Convert to a data.frame
    plt.table <- data.frame(fitted_values = fit_vals, pearson_residuals = pear_resid)

    # Ggplot
    fitted_vs_pearson <- base_plot +
      geom_point(
        data = plt.table,
        aes(x = .data[["fitted_values"]], y = .data[["pearson_residuals"]]),
        color = point.color, alpha = point.alpha, shape = point.shape
      ) +
      geom_hline(yintercept = ref.line.y.intercept, linetype = "dashed", color = "red", alpha = 0.5) +
      xlab("Fitted Values") + ylab("Pearson Residuals") +
      scale_x_continuous(breaks = round(seq(
        from = min(plt.table$fitted_values),
        to = max(plt.table$fitted_values), length.out = 10
      ))) +
      scale_y_continuous(breaks = round(seq(
        from = min(plt.table$pearson_residuals),
        to = max(plt.table$pearson_residuals), length.out = 6
      ))) +
      ggtitle("ScMaSigPro: Model Validation Plot", "Fitted Values vs Pearson Residual")

    # Plot
    return(fitted_vs_pearson)
  } else if (plot.type == "cooks_vs_obs") {
    # Ggplot
    p <- ggplot(glm.model, aes(seq_along(.cooksd), .cooksd)) +
      geom_bar(
        stat = "identity", position = "identity",
        width = bar_width, fill = point.color, alpha = point.alpha
      ) +
      xlab("Observations") +
      ylab("Cooks's Distance") +
      ggtitle("ScMaSigPro: Model Validation Plot", "Cook's Distance vs Observations") +
      theme_classic() +
      theme(
        panel.grid.major = element_line(size = 0.2, color = "lightgrey", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(0.7)),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
      ) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5)

    # Plot
    print(p)
  } else if (plot.type == "pearson_vs_obs") {
    # Extract pearson residuals
    pear_resid <- data.frame(
      obs = rownames(as.data.frame(resid(glm.model, type = "pearson"))),
      pearson_residual = resid(glm.model, type = "pearson")
    )

    if (func == "zeroinfl") {
      if (!is.null(input.data)) {
        obs_data <- input.data
        obs_data$observation <- rownames(obs_data)
      } else {
        stop("Supply Covriate Data")
      }
    } else if (func == "glm") {
      # Extraxt obsData fro each value
      obs_data <- as.data.frame(glm.model$data)

      # Rowname
      obs_data$observation <- rownames(obs_data)
    } else {
      stop("Function not supported")
    }

    # Melt
    obs_data <- melt(obs_data)

    # Rename Columns
    colnames(obs_data) <- c("obs", "covariate", "obs_value")

    # Merge to plot table
    plt.table <- merge(obs_data, pear_resid, "obs")

    # Ggplot
    pearson_vs_obs <- base_plot +
      geom_point(
        data = plt.table,
        aes(x = .data[["obs_value"]], y = .data[["pearson_residual"]]),
        color = point.color, alpha = point.alpha, shape = point.shape
      ) +
      geom_hline(yintercept = ref.line.y.intercept, linetype = "dashed", color = "red", alpha = 0.5) +
      xlab("Observed Values") + ylab("Pearson Residuals") +
      facet_wrap(~covariate)

    if (func == "zeroinfl") {
      pearson_vs_obs <- pearson_vs_obs + ggtitle("ScMaSigPro: Model Validation Plot", "Observed Value per covariate (Scaled) vs Pearson Residual")
    } else if (func == "glm") {
      pearson_vs_obs <- pearson_vs_obs + ggtitle("ScMaSigPro: Model Validation Plot", "Observed Value per covariate vs Pearson Residual")
    }

    # Plot
    return(pearson_vs_obs)
  } else if (plot.type == "fitted_vs_obs") {
    # Extract Table
    plt.table <- data.frame(
      obs_value = glm.model$y,
      fit_value = fitted(glm.model)
    )

    # Ggplot
    fitted_vs_obs <- base_plot +
      geom_point(
        data = plt.table, aes(x = .data[["fit_value"]], y = .data[["obs_value"]]),
        color = point.color, alpha = point.alpha, shape = point.shape
      ) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
      xlab("Fitted Values") + ylab("Observed Values") +
      scale_y_continuous(breaks = round(seq(
        from = 0,
        to = max(plt.table$obs_value),
        length.out = 10
      ))) +
      scale_x_continuous(breaks = round(seq(
        from = 0,
        to = max(plt.table$obs_value),
        length.out = 10
      )))

    ggtitle("ScMaSigPro: Model Validation Plot", "Fitted Values vs Observed Value")

    # Plot
    return(fitted_vs_obs)
  } else if (plot.type == "std_vs_fitted") {
    if (func == "zeroinfl") {
      # Convert to a data.frame
      plt.table <- data.frame(
        fitted_values = fitted(glm.model),
        std_residuals = resid(glm.model)
      )
    } else if (func == "glm") {
      # Convert to a data.frame
      plt.table <- data.frame(
        fitted_values = fitted(glm.model),
        std_residuals = rstandard(glm.model)
      )
    }

    # Add Layers
    std_vs_fitted <- base_plot +
      geom_point(
        data = plt.table,
        aes(x = .data[["fitted_values"]], y = .data[["std_residuals"]]),
        color = point.color, alpha = point.alpha, shape = point.shape
      ) +
      geom_hline(yintercept = ref.line.y.intercept, linetype = "dashed", color = "red", alpha = 0.5) +
      xlab("Fitted Values") + ylab("Standard Residuals") +
      scale_x_continuous(breaks = round(seq(
        from = min(plt.table$fitted_values),
        to = max(plt.table$fitted_values), length.out = 10
      ))) +
      scale_y_continuous(breaks = round(seq(
        from = min(plt.table$std_residuals),
        to = max(plt.table$std_residuals), length.out = 6
      ))) +
      ggtitle("ScMaSigPro: Model Validation Plot", "Fitted Values vs Standard Residuals")

    return(std_vs_fitted)
  } else {
    stop(paste(plot.type, "is not supported."))
  }
}
