#' sc_pVector
#'
#' @name sc_pVector
#' @aliases sc_pVector
#'
#' @param scMaSigPro.obj An object of the class scMaSigProClass should be supplied.
#' @param assay.name The name of the assay utilized from the Single Cell experiment class.
#' @param p.vector.sig The significance level of the p-value for full vs intercept models.
#' @param m.t.c.method  Method for multiple testing correction of p.value with  \code{\link[stats]{p.adjust}}
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
#' @param use.offset If set to true then log(Size factors) will be used as offsets
#'

#sc_pVector_test <- function(scMaSigPro.obj,
sc_pVector <- function(scMaSigPro.obj,
                            assay.name = "counts",
                            p.vector.sig = 0.05,
                       m.t.c.method = "BH",
                            min.counts = 6,
                            model.type = "nb",
                            model.control = list(ep = 0.00001, nb.k = 10, link = "log"),
                       use.offset = FALSE
                            ) {

    # Check Model.Type and Model.Control Consistancy
    assert_that(all(model.type %in% accpt.fam),
                msg = paste("'model.type' should be one of", accepted_model_string))
    assert_that(is.list(model.control),
               msg = "Invalid parameter: Expected a list for 'model.control'")

    # Set Seed
    set.seed(2022)

    # Prepare model control
    #source("model.control.R")
    model.parameters <- prepare.model.control(model_type = model.type,
                          model_control = model.control)

    # Covariate
    covariate <- as.data.frame(scMaSigPro.obj@covariate@covariate)

    # Extract Counts
    p.counts <- scMaSigPro.obj@assays@data@listData[[assay.name]]

    # Checking for NA
    p.counts <- p.counts[apply(p.counts, 1, count.na) >= min.counts, ]

    # Calculate Total sum
    # Drop genes if total sum is 0
    sum.tot <- apply(p.counts, 1, sum)
    p.counts.0 <- which(sum.tot == 0)
    if (length(p.counts.0) > 0) {
        p.counts <- p.counts[-p.counts.0, ]
    }
    if (use.offset == TRUE){
        size_factors <- estimateSizeFactorsForMatrix(p.counts+1)
        offset.data <- log(size_factors)
        print(sum(is.infinite(offset.data)))
        offset.data <- use.offset
    }

    # Start Parallel Compute
    #useCores <- makeCluster(detectCores(logical=FALSE)-1)
    #registerDoParallel(useCores)

    # Convert matrix to named list by row
    #computed_models <- parApply(useCores, p.counts, 1, function(gene.i.count) {
    computed_models <- apply(p.counts, 1, function(gene.i.count) {

        if (model.type %in% c("zinb", "zip")){

            if (min(gene.i.count) == 0){

            # Scale the Variables
            covariate <- apply(covariate, 2, scale)

            # Attach Gene expression to data
            covariate_data = cbind(gene.i.count, covariate)

            # Create Model formula Strings
            #"y~.|y~."
            y.full.zi.full.formula.str <- as.formula(paste(
                "gene.i.count", "~",
                paste(colnames(covariate_data)[!colnames(covariate_data) == "gene.i.count"],
                      collapse = "+"
                )))
            #"y~.|y~1"
            y.full.zi.intercept.formula.str <- as.formula(paste(
                "gene.i.count", "~",
                paste(colnames(covariate_data)[!colnames(covariate_data) == "gene.i.count"],
                      collapse = "+"
                ), "| 1"))

            #"y~1|y~."
            y.intercept.zi.full.formula.str <- as.formula(paste(
                "gene.i.count ~ 1 |",paste(
                    colnames(covariate_data)[!colnames(covariate_data) == "gene.i.count"],
                    collapse = "+")) )

            #"y~1|y~1"
            y.intercept.zi.intercept.formula.str <- as.formula("gene.i.count ~ 1 | 1")

            # Full Model #"y~.|y~."
            suppressWarnings(y.full.zi.full.model <- zeroinfl(
                formula = y.full.zi.full.formula.str,
                data = as.data.frame(covariate_data),
                dist = model.parameters$family
            ))

            # Intercept Model #"y~.|y~1"
            suppressWarnings(y.full.zi.intercept.model <- zeroinfl(
                formula = y.full.zi.intercept.formula.str,
                data = as.data.frame(covariate_data),
                dist = model.parameters$family
            ))

            # print(summary(y.full.zi.intercept.model))
            # stop("Expected Stop")

            # Y Inercept and Zi full # y.intercept.zi.full.formula.str
            suppressWarnings(y.intercept.zi.full.model <- zeroinfl(
                formula = y.intercept.zi.full.formula.str,
                data = as.data.frame(covariate_data),
                dist = model.parameters$family
            ))

            #print(summary(y.intercept.zi.full.model))
            #stop("Expected Stop")

            # Y Inercept and Zi full # y.intercept
            suppressWarnings(y.intercept.zi.intercept.model <- zeroinfl(
                formula = y.intercept.zi.intercept.formula.str,
                data = as.data.frame(covariate_data),
                dist = model.parameters$family
            ))

            # Create Model List
            tmp.model.list <- list(
                y_full_zi_full  = y.full.zi.full.model,
                y_intercept_zi_intercept  = y.intercept.zi.intercept.model,
                y_intercept_zi_full  = y.intercept.zi.full.model,
                y_full_zi_intercept  = y.full.zi.intercept.model
            )

            # Get the name
            selected.model.type <- getLowestAICModel(tmp.model.list)
            selected.model <- tmp.model.list[[selected.model.type]]

            # Test the selected model against the Intercept Model
            test.res <- pchisq(2 * (logLik(selected.model) - logLik(y.intercept.zi.intercept.model)),
                               df = 3, lower.tail = FALSE)
            zi.full.model.p.value <-  as.data.frame(test.res)[1,1]

            y.full.zi.full.model = selected.model
            #y.intercept.zi.full.model = y.intercept.zi.full.model

            }else{
                y.full.zi.full.model <- glm(gene.i.count ~ .,
                                  data = covariate#,
                                  #family = model.parameters$family,
                                  #epsilon = model.parameters$espsilon
                )

                # Compute intercept model
                y.intercept.zi.intercept.model <- glm(gene.i.count ~ 1,
                                       data = covariate#,
                                       #family = model.parameters$family,
                                       #epsilon = model.parameters$espsilon
                )

                # Test
                test.res <- anova(y.full.zi.full.model,
                                  y.intercept.zi.intercept.model,
                                  test = "Chisq"
                )
                # P value
                zi.full.model.p.value <- as.numeric(test.res[5][2, 1])
            }

            return(list(
                full.model = y.full.zi.full.model,
                intercept.model = y.intercept.zi.intercept.model,
                p.value = zi.full.model.p.value
            ))

        }else if (model.type %in% c("gnb", "gp")){

            # Attach Gene expression to data
            covariate_data = cbind(gene.i.count, covariate)

            # Create Model formula Strings
            full.formula.str <- as.formula(paste(
                "gene.i.count", "~",
                paste(colnames(covariate_data)[!colnames(covariate_data) == "gene.i.count"],
                      collapse = "+"
                )))
            intercept.formula.str <- as.formula("gene.i.count~1")

            # Full Model
            full.model <-glmmTMB(
                formula = full.formula.str,
                data = covariate_data,
                family = model.parameters$family
            )

            # Intercept Model
            intercept.model <-glmmTMB(
                formula = gene.i.count~1,
                data = covariate_data,
                family = model.parameters$family
            )

            # test of significance
            test.res <- anova(intercept.model,
                              full.model,
                              test = "Chisq"
            )

            full.model.p.value <- test.res[["Pr(>Chisq)"]][2]

            return(list(
                full.model = full.model,
                intercept.model = intercept.model,
                p.value = full.model.p.value
            ))
        }
        else if(model.type %in% stat.models){
            
            if (use.offset == TRUE){
                full.model <- glm(gene.i.count ~ .,
                                  data = covariate,
                                  family = model.parameters$family,
                                  epsilon = model.parameters$espsilon,
                                  offset = offset.data
                )
                
            }else if (use.offset == FALSE){
                # Compute full model
                full.model <- glm(gene.i.count ~ .,
                                  data = covariate,
                                  family = model.parameters$family,
                                  epsilon = model.parameters$espsilon
                )
            }

        # Compute intercept model
            if (use.offset == TRUE){
                intercept.model <- glm(gene.i.count ~ 1,
                                       data = covariate,
                                       family = model.parameters$family,
                                       epsilon = model.parameters$espsilon,
                                       offset = offset.data
                )
                
            }else if (use.offset == FALSE){
                intercept.model <- glm(gene.i.count ~ 1,
                                       data = covariate,
                                       family = model.parameters$family,
                                       epsilon = model.parameters$espsilon
                )
            }

        # MaSigPro, Pvector Evaluation
        test.res <- anova(intercept.model,
                          full.model,
                          test = "Chisq"
        )
        # P value
        full.model.p.value <- as.numeric(test.res[5][2, 1])

        # Return
        return(list(
            full.model = full.model,
            intercept.model = intercept.model,
            p.value = full.model.p.value
        ))

        }
    })

     # Stop Compute
    #stopCluster(cl = useCores)

    #--------------------------------------------------------------------------#
    #------------------ Extraction of Components for the Object----------------#
    #--------------------------------------------------------------------------#

    # Assign row names as names to each sublist
    names(computed_models) <- rownames(p.counts)

    # Extract Object
    full.model.list <- sapply(computed_models,  simplify = F, function(model.i){
        return( model.i[["full.model"]])
    })
    intercept.model.list <- sapply(computed_models, simplify = F, function(model.i){
        return(model.i[["intercept.model"]])
    })

    # Extract Pvalues
    p.value.vector <- sapply(computed_models, function(model.i) {
        return(model.i[["p.value"]])
    })

    # Adjusted P-Values
    adj.p.value.vector <- p.adjust(p.value.vector,
                                   method = m.t.c.method,
                                   n = length(p.value.vector)
    )

    # Update Parameters
    scMaSigPro.obj@parameters@assay.name = assay.name
    scMaSigPro.obj@parameters@p.vector.sig = p.vector.sig
    scMaSigPro.obj@parameters@m.t.c.method = m.t.c.method
    scMaSigPro.obj@parameters@min.counts = min.counts
    scMaSigPro.obj@parameters@model.type = model.type
    scMaSigPro.obj@parameters@model.control = model.control

    # Create Object
    pVector.obj <- new("pVectorClass",
                       full.models = full.model.list,
                       intercept.models = intercept.model.list,
                       p.value = p.value.vector,
                       adj.p.value = adj.p.value.vector
    )

    #--------------------------------------------------------------------------#
    #------------------ Adding and returning the objects ----------------------#
    #--------------------------------------------------------------------------#

    # Add to object
    scMaSigPro.obj@pVector <- pVector.obj

    # Return adjusted p-value vector
    return(scMaSigPro.obj)
}
