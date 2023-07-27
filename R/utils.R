# Count na
count.na <- function(x) {
  return(length(x) - length(x[is.na(x)]))
}

# Accepted families
accpt.fam <- c("g", "nb", "zinb", "p", "zip", "zap", "zanb", "gp", "gnb")
stat.models <- c("g", "nb", "p")
glmmtmb.models <- c("zinb","zip", "zap", "zanb", "gp", "gnb")

# String Messages
accepted_model_string <- '"g": Gaussian, "nb",
"zinb": Zero-Inflated Negative Binomial,
"p": Poisson, "zip": Zero-Inflated Poisson,
"zap": Zero-Altered Poisson,
"zanb": Zero-Altered Negative Binomial,
"gp": Generalized Poisson,
"gnb": Generalized Negative Binomial'

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

create_range <- function(x) {
    # Convert the factor into Character
    y <- as.character(x[["bin"]])
    y <- y %>% str_remove_all(pattern = "\\[|\\]|\\(|\\)")
    y1 <- as.numeric(sapply(strsplit(y, ","), "[", 1))
    y2 <- as.numeric(sapply(strsplit(y, ","), "[", 2))
    rangeVec <- c(y1, y2, x[["bin_size"]], x[["binned_time"]])
    return(as.numeric(rangeVec))
}

calc_bin_size <- function(x) {
    size <- length(
        c(stringr::str_split(x[["cluster.members"]], "\\|"))[[1]]
    )
    return(as.numeric(size))
}
