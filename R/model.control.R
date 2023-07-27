# Need a check function for all the swicth statements
prepare.model.control <- function(model_type, model_control) {
  result <- switch(model_type,
    "nb" = {
      fam.obj <- negative.binomial(
        theta = model_control$nb.k,
        link = model_control$link
      )
      espsilon <- model_control$ep
      return(list(
        family = fam.obj, espsilon = espsilon,
        theta = model_control$nb.k,
        link = model_control$link
      ))
    },
    "p" = {
      fam.obj <- poisson(link = model_control$link)
      return(list(
        family = fam.obj,
        link = model_control$link
      ))
    },
    "g" = {
      fam.obj <- gaussian(link = model_control$link)
      return(list(
        family = fam.obj,
        link = model_control$link
      ))
    },
    "gp" = {
      fam.obj <- genpois(link = model_control$link)
      return(list(
        family = fam.obj,
        link = model_control$link
      ))
    },
    "gnb" = {
      fam.obj <- nbinom2(link = model_control$link)
      return(list(
        family = fam.obj,
        link = model_control$link
      ))
    },
    "zip" = {
      fam.obj <- "poisson"
      return(list(
        family = fam.obj,
        link = NULL
      ))
    },
    "zinb" = {
      fam.obj <- "negbin"
      return(list(
        family = fam.obj,
        link = NULL
      ))
    },
    "default" = {
      print("Supplied Negative Binomial")
    }
  )
}
