#' @title Plot the diagnostics of the model of a gene
#'
#' @description
#' `plotDiagnostics()` draws four diagnostic plots of a model of a gene of
#' the gene lists of [sc.filter()]: the normal Q-Q plot of the standardised
#' residuals, the standardised residuals and their square root against the
#' fitted values (scale-location), and the observed against the predicted
#' counts.
#'
#' @param model The model: `"optimized"`, the model selected by [sc.t.fit()];
#'   `"full"`; or `"intercept"`, the intercept-only model.
#'
#' @importFrom stats residuals fitted rstandard predict
#' @importFrom ggplot2 .data aes element_blank element_line geom_abline
#'   geom_hline geom_point geom_smooth ggplot labs stat_qq stat_qq_line theme
#'   theme_minimal
#'
#' @return A patchwork object with four ggplot panels.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @family plots
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'     plotDiagnostics(scmp.ob, "Gene10")
#' }
#'
#' @export
plotDiagnostics <- function(scmpObj,
                            feature_id,
                            model = "optimized") {
    # Check Object Validity
    .check_scmp(scmpObj)

    # The fits of the current design: a step that changes the bins, the
    # counts or the design empties them, and an object saved by an earlier
    # version after sc.set.poly() keeps fits without their offset
    fitted <- length(scmpObj@Profile@p_values) > 0L &&
        (!scmpObj@Parameters@offset || length(scmpObj@Design@offset) > 0L)
    .check(fitted,
        msg = "The current design has not been fitted: run sc.p.vector() first."
    )

    # Check if the results exist
    .check(!isEmpty(scmpObj@Significant@genes),
        msg = paste("No Significant genes found, please run the workflow first")
    )

    # Check model type
    model <- .check_choice(model, c("optimized", "full", "intercept"),
        msg = paste(
            "Invalid 'model', please select one of 'optimized', 'full' or",
            "'intercept'"
        )
    )

    # The four panels are arranged with the operators of patchwork
    .check(requireNamespace("patchwork", quietly = TRUE),
        msg = paste(
            "plotDiagnostics() needs the 'patchwork' package to arrange its",
            "panels; install it with install.packages('patchwork')."
        )
    )

    # Get all the available genes
    avail_genes <- unique(unlist(scmpObj@Significant@genes))

    # Check if the requested gene is availble as singificant gene
    .check(all(feature_id %in% avail_genes),
        msg = paste(
            "The requested gene is not available in the significant genes",
            "list"
        )
    )

    # The counts of the gene in the bins of the design
    dis <- scmpObj@Design@predictor_matrix
    .check_dense_counts(scmpObj, bins = rownames(dis))
    y_mat <- as.matrix(
        assay(scmpObj@Dense, "bulk.counts", withDimnames = FALSE)
    )
    y <- y_mat[feature_id, rownames(dis)]

    # Set theme for the plot
    theme_set <- theme_minimal(base_size = 12) +
        theme(
            legend.position = "bottom",
            panel.grid.major = element_line(
                color = "grey90", linewidth = 0.3, linetype = "dashed"
            ),
            panel.grid.minor = element_blank()
        )

    # The terms of the model: "full" and "intercept" use the design,
    # "optimized" the terms that sc.t.fit() selected for the gene, those
    # with a p-value in showSol() (its columns after the first three); a
    # model without terms is the intercept-only model
    if (model == "optimized") {
        selected <- !is.na(unlist(showSol(scmpObj)[feature_id, -seq_len(3)]))
        if (any(selected)) {
            dis <- dis[, selected, drop = FALSE]
        } else {
            model <- "intercept"
        }
    }

    # The counts and the terms, named y and x1, x2, ... by position as in
    # the fits (see .fit_frame()), so that any design names work
    model_df <- data.frame(y = y, .fit_frame(dis))

    # Extract offsets
    if (scmpObj@Parameters@offset) {
        offset_data <- scmpObj@Design@offset
    } else {
        offset_data <- NULL
    }

    # The gene's own theta when sc.p.vector() estimated one per gene
    family <- .gene_family(scmpObj@Parameters@distribution, feature_id)

    if (model == "intercept") {
        glm_s3 <- glm(y ~ 1,
            data = model_df,
            family = family,
            epsilon = scmpObj@Parameters@epsilon,
            offset = offset_data,
            weights = NULL,
            maxit = scmpObj@Parameters@max_it
        )
    } else {
        # Fit Model
        glm_s3 <- glm(y ~ .,
            data = model_df,
            family = family,
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
