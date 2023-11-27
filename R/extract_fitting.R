#' @title Extract Fitting
#' 
#' @param reg 
#' @param lmf 
#' @param model.glm.0 
#' @param dis 
#' @param family
#' @param name 
#' @param vars.in 
#' @param alfa 
#' @param influ.info
#' 
#' @keywords internal

extract_fitting <- function(reg, lmf, model.glm.0, dis, family, name, vars.in, alfa, influ.info) {
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  y <- reg$y
  result <- summary(lmf)
  novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 4]))]
  influ <- influence.measures(reg)$is.inf
  influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
  influ1 <- which(apply(influ, 1, all))
  if (length(influ1) != 0) {
    paste.names <- function(a) {
      paste(names(a)[a], collapse = "/")
    }
    match <- match(rownames(dis), rownames(influ))
    influ <- as.data.frame(apply(influ, 1, paste.names))
    influ.info <- cbind(influ.info, influ[match, ])
    colnames(influ.info)[ncol(influ.info)] <- name
    influ.info <- as.matrix(influ.info)
  }
  result <- summary(reg)
  if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == "gaussian") | family$family != "gaussian") {
    # k <- i

    # Computing p-values
    # model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon, offset = offsetData)

    if (family$family == "gaussian") {
      test <- anova(model.glm.0, reg, test = "F")
      p.value <- test[6][2, 1]
    } else {
      test <- anova(model.glm.0, reg, test = "Chisq")
      p.value <- test[5][2, 1]
    }
    # Computing goodness of fitting:

    bondad <- (reg$null.deviance - reg$deviance) / reg$null.deviance
    if (bondad < 0) {
      bondad <- 0
    }
    beta.coeff <- result$coefficients[, 1]
    beta.p.valor <- result$coefficients[, 4]
    coeff <- rep(0, (length(vars.in) + 1))
    if (length(novar) != 0) {
      for (m in 1:length(novar)) {
        coeff[position(dis, novar[m]) + 1] <- NA
      }
    }
    p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 1)))

    if (result$coefficients[, 4][rownames(result$coefficients) ==
      "(Intercept)"] < alfa) {
      coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) ==
        "(Intercept)"]
      p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) ==
        "(Intercept)"]
      t[1] <- result$coefficients[, 3][rownames(result$coefficients) ==
        "(Intercept)"]
    }
    for (j in 2:length(coeff)) {
      if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
        coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) ==
          vars.in[j - 1]]
        p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) ==
          vars.in[j - 1]]
        t[j] <- result$coefficients[, 3][rownames(result$coefficients) ==
          vars.in[j - 1]]
      }
    }
    if (!all(is.na(p.valor))) {
      sol <- rbind(sol, as.numeric(c(
        p.value, bondad,
        p.valor
      )))
      coefficients <- rbind(coefficients, coeff)
      t.score <- rbind(t.score, t)
      sig.profiles <- rbind(sig.profiles, y)
      h <- nrow(sol)
      rownames(sol)[h] <- name
      rownames(coefficients)[h] <- name
      rownames(t.score)[h] <- name
      rownames(sig.profiles)[h] <- name
    }
  }

  # Return Calculation
  return(list(
    p_value = p.value,
    bondad = bondad,
    p_valor = p.valor,
    coeff = coeff,
    t = t,
    sig_profiles = y,
    sol = sol,
    influ.info = influ.info,
    feature_name = name
  ))
}
