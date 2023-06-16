#' sc_tFit
#'
#' @name sc_tFit
#' @aliases sc_tFit
#'
#' @param scMaSigPro.obj An object of the class scMaSigProClass should be supplied.
#' @param covar.sig The significance level of the p-value for each covariate during model selection.
#' @param min.counts Genes with less than this number of true numerical values will be excluded from the analysis.
#' @param model.type One of the  "g": Gaussian, "nb", "zinb": Zero-Inflated Negative Binomial,
#' "p": Poisson, "zip": Zero-Inflated Poisson, "zap": Zero-Altered Poisson, "zanb": Zero-Altered Negative Binomial,
#' "gp": Generalized Poisson, "gnb": Generalized Negative Binomial'.
#' @param model.control A list of control parameters for the model estimation. For e.g.
#'  The available control parameters for 'nb': Negative Binomial are:
#' \describe{
#'      \item{\code{ep}}{The epsilon value for model convergence. Default is 0.00001.}
#'      \item{\code{nb.k}}{The size parameter for the negative binomial distribution. Default is 10.}
#'   }
#'

sc_tFit <- function(scMaSigPro.obj,
                    covar.sig = 0.05,
                    min.counts = 6,
                    model.type = "nb",
                    model.control = list(ep = 0.00001, nb.k = 10, link = "log")
) {

    # Check Model.Type and Model.Control Consistancy
    assert_that(all(model.type %in% accpt.fam),
                msg = paste("'model.type' should be one of", accepted_model_string))
    assert_that(is.list(model.control),
                msg = "Invalid parameter: Expected a list for 'model.control'")

    # Set Seed
    set.seed(2022)

    # Make Family
    model.parameters <- prepare.model.control(model_type = model.type,
                                              model_control = model.control)

  # Identify significant genes based on adjusted p-value
  sig.genes <- names(scMaSigPro.obj@pVector@adj.p.value[scMaSigPro.obj@pVector@adj.p.value <= scMaSigPro.obj@parameters@p.vector.sig])

  # Get models associated with significant genes
  sig.full.models.list <- scMaSigPro.obj@pVector@full.models[sig.genes]
  sig.intercept.models.list <- scMaSigPro.obj@pVector@intercept.models[sig.genes]

  # Extract the models that are significant
  sig.models.list <- mapply(
    function(x, y) {
      list("full.models" = x, "intercept.models" = y)
    },
    sig.full.models.list, sig.intercept.models.list,
    SIMPLIFY = FALSE
  )

  # Compute for STAT Models
  if (scMaSigPro.obj@parameters@model.type %in% c("g", "p", "nb")){

  # Apply a function to each significant gene model
  sig.covars.list <- lapply(sig.models.list, function(gene.model.i, USE.NAMES = TRUE, dist_fam = dist, epsilon = ep) {
    # Extract the Intercept Model and Full Model
    intercept.model <- gene.model.i[["intercept.models"]]
    covar.model <- gene.model.i[["full.models"]]

    # Obtain a full model summary
    covar.model.summary <- summary(covar.model)

    # Get the maximum p-value out of the full model
    # The pvalue could be NA
    covar.model.covariate.max.pvalue <- max(covar.model.summary$coefficients[, 4][-1], na.rm = TRUE)

    # Extract the model frame for iterative model fitting
    covar.model.data <- model.frame(covar.model)

    # Check if the any of the pvalue is greater than the asked significance level
    if (covar.model.covariate.max.pvalue > covar.sig) {
      # Initiate a while-loop
      while (covar.model.covariate.max.pvalue > covar.sig) {
        # Identify the covariate with the maximum p-value
        coefficient.table <- data.frame(covariate_p_value = covar.model.summary$coefficients[, 4][-1])

        # Get the Irrelevant.covariate
        # Check Required for NA and Null
        irrelevant.covariate <- rownames(coefficient.table[coefficient.table[["covariate_p_value"]] == covar.model.covariate.max.pvalue, , drop = F])

        # Remove the identified covariate from the coefficient.table
        covar.model.data <- covar.model.data[, !(colnames(covar.model.data) == irrelevant.covariate), drop = F]

        # Check if there are still covariates in the model
        if (ncol(covar.model.data) >= 2) {
          # Create Formula
          formula.str <- paste(
            "gene.i.count", "~",
            paste(colnames(covar.model.data)[!colnames(covar.model.data) == "gene.i.count"],
              collapse = "+"
            )
          )

          # Refit the model without the removed covariate
          covar.model <- glm(
            formula = formula.str,
            data = covar.model.data,
            family = model.parameters$family,
            epsilon = model.parameters$espsilon
          )

          # Extract summary
          covar.model.summary <- summary(covar.model)

          # Extract the Coefficient table
          coefficient.table <- data.frame(covariate_p_value = covar.model.summary$coefficients[, 4][-1])

          # Extract the max-pvalue
          # Check Required could be NA or NULL
          covar.model.covariate.max.pvalue <- max(coefficient.table[["covariate_p_value"]], na.rm = TRUE)
        } else if (ncol(covar.model.data) <= 1) {
          # Transfer the Intercept Model
          covar.model <- intercept.model

          # Exit with Zero
          covar.model.covariate.max.pvalue <- 0
        }
      }

      # Return the Selected Model
      return(covar.model)
    } else if (covar.model.covariate.max.pvalue <= covar.sig) {
      return(covar.model)
    }
  })
  }else if (scMaSigPro.obj@parameters@model.type %in% c("zip", "zinb")){


      sig.covars.list <- lapply(sig.models.list, function(gene.model.i, USE.NAMES = TRUE, dist_fam = model.control$family) {

          intercept.model <- gene.model.i[["intercept.models"]]
          full.model <- gene.model.i[["full.models"]]
          covar.model <- full.model

          # Extract Summary of the full models summary and Count
          full.model.summary <- summary(full.model)

          # Fixing the Count Component
          count.full.model.coeff <- as.data.frame(full.model.summary[["coefficients"]][["count"]])

          # Get the maximum p-value out of the full model
          # The pvalue could be NA
          covar.model.covariate.max.pvalue <- max(count.full.model.coeff[, 4][-1], na.rm = TRUE)

          # Extract the model frame for iterative model fitting
          covar.model.data <- model.frame(full.model)

          # Check if the any of the pvalue is greater than the asked significance level
          if (covar.model.covariate.max.pvalue > covar.sig) {

              # Initiate a while-loop
              while (covar.model.covariate.max.pvalue > covar.sig) {
                  # Identify the covariate with the maximum p-value
                  coefficient.table <- as.data.frame(full.model.summary[["coefficients"]][["count"]])
                  colnames(coefficient.table) <- c("Estimate", "Std_Err", "Z_Value", "covariate_p_value")

                  # Get the Irrelevant.covariate
                  # Check Required for NA and Null
                  irrelevant.covariate <- rownames(coefficient.table[coefficient.table[["covariate_p_value"]] == covar.model.covariate.max.pvalue, , drop = F])

                  # Remove the identified covariate from the coefficient.table
                  covar.model.data <- covar.model.data[, !(colnames(covar.model.data) == irrelevant.covariate), drop = F]

                  # Check if there are still covariates in the model
                  if (ncol(covar.model.data) >= 2) {
                      # Create Formula
                      formula.str <- as.formula(paste(
                          "gene.i.count", "~", paste(colnames(covar.model.data)[!colnames(covar.model.data) == "gene.i.count"],
                                                     collapse = "+"
                                                     )
                          ))

                      # Refit the model without the removed covariate
                      suppressWarnings(covar.model <- zeroinfl(
                          formula = formula.str,
                          data = as.data.frame(covar.model.data),
                          dist = dist_fam
                      ))

                      # Extract summary
                      covar.model.summary <- summary(covar.model)

                      # Extract the Coefficient table
                      coefficient.table <- data.frame(covar.model.summary[["coefficients"]][["count"]])

                      # Extract the max-pvalue
                      # Check Required could be NA or NULL
                      covar.model.covariate.max.pvalue <- max(coefficient.table[, 4][-1], na.rm = TRUE)

                      } else if (ncol(covar.model.data) <= 1) {
                      # Transfer the Intercept Model
                      covar.model <- intercept.model

                      # Exit with Zero
                      covar.model.covariate.max.pvalue <- 0
                  }
              }

              # Fix the Zero Inflation

              # Extract Summary of the full models summary and Count
              covar.model.summary <- summary(covar.model)

              # Fixing the Count Component
              zero.full.model.coeff <- as.data.frame(covar.model.summary[["coefficients"]][["zero"]])

              # Get the maximum p-value out of the full model
              # The pvalue could be NA
              covar.model.covariate.max.pvalue <- max(zero.full.model.coeff[, 4][-1], na.rm = TRUE)

              # Extract the model frame for iterative model fitting
              covar.model.data <- model.frame(covar.model)
              drop.covars <- colnames(covar.model.data)[colnames(covar.model.data) != "gene.i.count"]

              # Check if the any of the pvalue is greater than the asked significance level
              if (covar.model.covariate.max.pvalue > covar.sig) {

                  # Initiate a while-loop
                  while (covar.model.covariate.max.pvalue > covar.sig) {

                      # Identify the covariate with the maximum p-value
                      coefficient.table <- as.data.frame(covar.model.summary[["coefficients"]][["zero"]])
                      colnames(coefficient.table) <- c("Estimate", "Std_Err", "Z_Value", "covariate_p_value")

                      # Get the Irrelevant.covariate
                      # Check Required for NA and Null
                      irrelevant.covariate <- rownames(coefficient.table[coefficient.table[["covariate_p_value"]] == covar.model.covariate.max.pvalue, , drop = F])
                      drop.covars <- drop.covars[!(drop.covars %in% irrelevant.covariate)]

                      # Remove the identified covariate from the coefficient.table
                      #covar.model.data <- covar.model.data[, !(colnames(covar.model.data) == irrelevant.covariate), drop = F]

                      if (length(drop.covars) > 0){
                      # Check if there are still covariates in the model
                          # Create Formula
                          formula.str <- as.formula(paste(
                              "gene.i.count ~ .|", paste(drop.covars, collapse = "+" )
                          ))

                          # Refit the model without the removed covariate
                          suppressWarnings(covar.model <- zeroinfl(
                              formula = formula.str,
                              data = as.data.frame(covar.model.data),
                              dist = dist_fam
                          ))

                          # Extract summary
                          covar.model.summary <- summary(covar.model)

                          # Extract the Coefficient table
                          coefficient.table <- data.frame(covar.model.summary[["coefficients"]][["zero"]])

                          # Extract the max-pvalue
                          # Check Required could be NA or NULL
                          covar.model.covariate.max.pvalue <- max(coefficient.table[, 4][-1], na.rm = TRUE)

                          if (length(drop.covars) == 1) {
                              covar.model <- covar.model
                              covar.model.covariate.max.pvalue <- 0
                          }
                      }
                  }
              }

              # Return the Selected Model
              return(covar.model)
          } else if (covar.model.covariate.max.pvalue <= covar.sig) {
              return(covar.model)
          }
      })

  }

  # Extract the coefficents and P-values
  measure_list <- lapply(names(sig.models.list), function(gene.model.i, USE.NAMES = TRUE, dist_fam =  scMaSigPro.obj@parameters@model.type, epsilon = ep,
                                                          sig_models_list = scMaSigPro.obj@pVector@intercept.models,
                                                          sig_covars_list = sig.covars.list,
                                                          covariate.vector = scMaSigPro.obj@covariate@covariate.vector,
                                                          group.vector = scMaSigPro.obj@covariate@group.vector) {

    # Extract Covar Model
    covar.model <- sig_covars_list[[gene.model.i]]

    # Extract the Intercept model
    intercept.model <- sig_models_list[[gene.model.i]]

    # For Stat Models
    if (dist_fam %in% c("g", "p", "nb")){

        # Calculate the test
        annova.chisq <- annova.chisq <- anova(intercept.model,
                                              covar.model,
                                              test = "Chisq"
        )

        # Extract the overall model p-value
        model.p.value <- as.numeric(annova.chisq[["Pr(>Chi)"]][2])

        # Compute Variance/ Goodness of fit
        model.r.square <- as.numeric((covar.model$null.deviance - covar.model$deviance) / covar.model$null.deviance)


        if (model.r.square < 0) {
            model.r.square <- 0
        }

        # Extract Summary as data frame
        covar.model.frame <- data.frame(summary(covar.model)$coefficients)

    }else if (dist_fam %in% c("zip", "zinb")){

        # Calculate Model pvalue
        test.res <- pchisq(2 * (logLik(covar.model) - logLik(intercept.model)), df = 3, lower.tail = FALSE)
        model.p.value <-  as.data.frame(test.res)[1,1]

        # Compute log-likelihoods
        LL_null <- logLik(intercept.model)
        LL_full <- logLik(covar.model)

        # Compute the Cox & Snell R^2
        R2_CS <- 1 - exp((2/nrow(covar.model$model)) * (LL_null - LL_full))

        # Compute the Nagelkerke's R^2
        R2_Nagelkerke <- R2_CS / (1 - (exp((2 * LL_null) / nrow(covar.model$model))))

        # Extract Pseudo R2
        model.r.square <- R2_Nagelkerke[1]

        # Extract Summary as data frame
        covar.model.frame <- data.frame(summary(covar.model)$coefficients$count)
    }

    # Get the significant Covariates
    nonsignifi.covariates <- c(
      setdiff(covariate.vector, rownames(covar.model.frame)[-1])
    )

    # Check the
    if (length(nonsignifi.covariates) > 0) {
      # Create an empty dataframe of nonsignifi.covariates
      nonsignifi.covariates.df <- data.frame(matrix(NA,
        nrow = length(nonsignifi.covariates),
        ncol = ncol(covar.model.frame)
      ))
      colnames(nonsignifi.covariates.df) <- colnames(covar.model.frame)
      rownames(nonsignifi.covariates.df) <- nonsignifi.covariates
      covar.model.frame <- rbind(covar.model.frame, nonsignifi.covariates.df)
      covar.model.frame <- covar.model.frame[c("(Intercept)", covariate.vector), , drop = F]
    }

    # Rename dimensions
    rownames(covar.model.frame) <- c("beta_0", paste0("beta_", rownames(covar.model.frame)[-1]))

    if (dist_fam %in% c("p", "nb", "g")){colnames(covar.model.frame) <- c("Estimate", "Std_Err", "t_value", "p_value")
    }else if(dist_fam %in% c("zip", "zinb")){
        colnames(covar.model.frame) <- c("Estimate", "Std_Err", "z_value", "p_value")
    }
    covar.model.frame[is.na(covar.model.frame[["Estimate"]]), "Estimate"] <- 0

    # Transpose the Dataframe
    covar.model.frame <- t(covar.model.frame)

    # Create named vector for p_values
    sol <- c(model.p.value, model.r.square, c(covar.model.frame["p_value", , drop = T]))
    names(sol) <- c("model.pvalue", "model.rsquare", colnames(covar.model.frame))

    # Create named vector for estimates
    beta_estimates <- c(covar.model.frame["Estimate", , drop = T])
    names(beta_estimates) <- colnames(covar.model.frame)

    # Create named vector for t_scores
    if (dist_fam %in% c("p", "nb", "g")){
        t_value <- c(covar.model.frame["t_value", , drop = T])
    }else if(dist_fam %in% c("zip", "zinb")){
        t_value <- c(covar.model.frame["z_value", , drop = T])
    }
    names(t_value) <- colnames(covar.model.frame)

    # Return
    return(list(
      sol = sol,
      beta_estimates = beta_estimates,
      t_value = t_value
    ))
  })

  names(measure_list) <- names(sig.models.list)

  # Update Parameters
  scMaSigPro.obj@parameters@covar.sig = covar.sig

  # Create Tfit Object
  tFit.obj <- new("tFitClass",
    covar.models = sig.covars.list,
    model.attributes = measure_list
  )

  # Add to MaSigPro Object
  scMaSigPro.obj@tFit <- tFit.obj

  return(scMaSigPro.obj)
}
