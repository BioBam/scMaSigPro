#' @keywords internal
# Accepted families
accpt.fam <- c("g", "nb", "zinb", "p", "zip", "zap", "zanb", "gp", "gnb")

#' @keywords internal
# Statistical models
stat.models <- c("g", "nb", "p")

#' @keywords internal
# Generalized linear mixed models
glmmtmb.models <- c("zinb", "zip", "zap", "zanb", "gp", "gnb")

#' @keywords internal
# String Messages
accepted_model_string <- '"g": Gaussian, "nb",
"zinb": Zero-Inflated Negative Binomial,
"p": Poisson, "zip": Zero-Inflated Poisson,
"zap": Zero-Altered Poisson,
"zanb": Zero-Altered Negative Binomial,
"gp": Generalized Poisson,
"gnb": Generalized Negative Binomial'
