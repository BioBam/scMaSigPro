# ScMaSigProClass with all the super-classes
# Author: Priyansh Srivastava
# Date: 2019-05-01

###############################################################################
#' @title ParameterConfig
#'
#' @description
#' S4 Class to hold parameters used in the analysis. It is constantly updated
#' asthe analysis progresses.
#'
#' @slot bin_ptime_col A character string representing the column name
#' for binned pseudotime values in dense data.
#' @slot path_prefix A character string representing the prefix for path labeling.
#' @slot root_label A character string representing the label for the tree root.
#' @slot pseudotime_colname A character string representing the column name for pseudotime values.
#' @slot bin_method A character string representing the algorithm used for binning.
#' @slot path_colname A character string representing the column name for path values.
#' @slot bin_colname A character string representing the column name for bin values.
#' @slot bin_size_colname A character string representing the column name for bin sizes.
#' @slot bin_members_colname A character string representing the column name for bin members.
#' @slot annotation_col A character string representing the column name for cell type annotation.
#' @slot g Integer. Number of genes taken in the regression fit.
#' @slot Q Numeric. Significance Level.
#' @slot min.na Numeric. Minimum value to estimate the model (degree+1) x Groups + 1. (Default = 6).
#' @slot MT.adjust A character string specifying the Pvalue correction method used.
#' @slot epsilon Numeric. Convergence tolerance.
#' @slot step.method A character string specifying the imputed step method for stepwise regression.
#' @slot offset A logical value specifying whether to use offset during model fitting.
#' @slot logOffset A logical value specifying whether to take the logarithm of the offsets during model fitting.
#' @slot max_it Integer. Maximum number of iterations to fit the model.
#' @slot poly_degree Integer with the polynomial degree to fit the regression. 1
#' @slot distribution Distribution used
#' @slot cluster.method Description
#' @slot cluster.fill.dim description
#' @slot cluster.fill.na description
#' specifies a linear regression, 2 a quadratic regression, etc.
#'
#' @name ParameterConfig
#' @aliases ParameterConfig-class
#' @rdname ParameterConfig-class
#' @importFrom methods is new
#' @keywords classes

setClass(
  "ParameterConfig",
  representation(
    bin_ptime_col = "character",
    path_prefix = "character",
    root_label = "character",
    pseudotime_colname = "character",
    step.method = "character",
    bin_method = "character",
    path_colname = "character",
    bin_colname = "character",
    bin_size_colname = "character",
    bin_members_colname = "character",
    annotation_col = "character",
    g = "integer",
    Q = "numeric",
    min.na = "numeric",
    MT.adjust = "character",
    epsilon = "numeric",
    offset = "logical",
    logOffset = "logical",
    max_it = "integer",
    poly_degree = "integer",
    distribution = "ANY",
    cluster.method = "character",
    cluster.fill.dim = "character",
    cluster.fill.na = "character"
  ),
  validity = function(object) {
    errors <- character(0)

    # Check if any character slots are empty or not of type character
    char_slots <- c(
      "bin_ptime_col", "path_prefix", "root_label",
      "pseudotime_colname", "bin_method",
      "path_colname", "bin_colname", "bin_size_colname",
      "bin_members_colname", "MT.adjust", "step.method",
      "annotation_col", "cluster.method", "cluster.fill.dim",
      "cluster.fill.na"
    )

    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name) # Corrected line
      if (length(slot_value) == 0 || !is.character(slot_value)) { # Corrected line
        errors <- c(errors, paste(slot_name, "must not be empty and should be of type character."))
      }
    }

    # Check for slot g
    if (!is.integer(object@g)) {
      stop("Slot 'g' must be an integer.")
    }


    # Check for slot Q
    if (!is.numeric(object@Q)) {
      stop("Slot 'Q' must be numeric.")
    }

    if (!is.numeric(object@epsilon)) {
      stop("Slot 'epsilon' must be numeric.")
    }

    # Check for slot min.na
    if (!is.numeric(object@min.na)) {
      stop("Slot 'min.na' must be an integer.")
    }

    # Check if any of the character slots have multiple values
    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name) # Corrected line
      if (length(slot_value) > 1) { # Corrected line
        errors <- c(errors, paste(slot_name, "should not contain multiple values."))
      }
    }

    if (length(errors) == 0) TRUE else errors
  },
  prototype = list(
    bin_ptime_col = "scmp_binned_pseudotime",
    path_prefix = "Path",
    root_label = "root",
    pseudotime_colname = "Pseudotime",
    path_colname = "Path",
    bin_method = "Sturges",
    bin_colname = "scmp_bin",
    g = 0L,
    Q = 0.05,
    min.na = 6,
    bin_size_colname = "scmp_bin_size",
    bin_members_colname = "scmp_bin_members",
    MT.adjust = "BH",
    step.method = "backward",
    epsilon = 1e-8,
    offset = TRUE,
    logOffset = FALSE,
    max_it = 100L,
    annotation_col = "cell_type",
    poly_degree = 2L,
    distribution = MASS::negative.binomial(10),
    cluster.method = "hclust",
    cluster.fill.dim = "col",
    cluster.fill.na = "zero"
  )
)

###############################################################################
#' @title MatrixDesign
#'
#' @description
#' S4 class to store Branching Path Assignments and Prediction Matrix
#'
#' @slot predictor A data frame containing the design matrix for the experiment.
#' @slot groups.vector A character vector specifying the experimental group to which
#' each variable belongs.
#' @slot alloc A data frame describing the detailed experimental design. Rows must contain
#' cells and columns experiment descriptors. The matrix must be binarized.
#'
#' @name MatrixDesign
#' @aliases MatrixDesign-class
#' @rdname MatrixDesign-class
#' @importFrom methods is new
#' @keywords classes

#'
#' @section Validity:
#'   Valid objects must have:
#'   \itemize{
#'     \item{predictor}{A valid data frame.}
#'     \item{groups.vector}{A valid character vector.}
#'     \item{alloc}{A valid data frame.}
#'   }
#'

setClass(
  "MatrixDesign",
  representation(
    predictor = "matrix",
    groups.vector = "character",
    alloc = "matrix"
  ),
  prototype = list(
    predictor = matrix(NA, nrow = 0, ncol = 0),
    groups.vector = character(),
    alloc = matrix(NA, nrow = 0, ncol = 0)
  ),
  validity = function(object) {
    if (!validObject(object@predictor)) {
      stop("predictor slot is not a valid data frame.")
    }
    if (!validObject(object@groups.vector)) {
      stop("groups.vector slot is not a valid character vector.")
    }
    if (!validObject(object@alloc)) {
      stop("poly_degree slot is not a valid matrix")
    }
    TRUE
  }
)
###############################################################################
#' @title VariableProfiles
#'
#' @description
#' S4 class to store results of the global model fitting for all the features.
#'
#' @slot non.flat a charcter vector
#' @slot p.vector Numeric vector containing the computed p-values.
#' @slot p.adjusted Numeric vector of FDR-adjusted p-values.
#' @slot FDR P-value at FDR \code{Q} control when Benjamini & Hochberg (BH) correction is used.
#'
#' @name VariableProfiles
#' @aliases VariableProfiles-class
#' @rdname VariableProfiles-class
#' @keywords classes

#' @author Priyansh Srivastava <spriyansh29@@gmail.com>
#' @seealso \code{\link{T.fit}}, \code{\link{lm}}
#' @importFrom stats family gaussian poisson
#' @importFrom utils data combn
#' @importFrom MASS negative.binomial

# Define the VariableProfiles with the following slots:
setClass("VariableProfiles",
  slots = c(
    non.flat = "character", # Matrix containing the expression values for significant genes
    p.vector = "numeric", # Matrix containing the computed p-values
    p.adjusted = "numeric", # Vector of FDR-adjusted p-values
    FDR = "numeric" # P-value at FDR Q control when Benjamini & Hochberg (BH) correction is used
  ),
  validity = function(object) {
    # Check for slot SELEC
    if (!is.character(object@non.flat)) {
      stop("Slot 'non.flat' must be a character")
    }

    # Check for slot sc.p.vector
    if (!is.numeric(object@p.vector)) {
      stop("Slot 'sc.p.vector' must be a numeric")
    }

    # Check for slot p.adjusted
    if (!is.numeric(object@p.adjusted)) {
      stop("Slot 'p.adjusted' must be numeric.")
    }

    # Check for slot FDR
    if (!is.numeric(object@FDR)) {
      stop("Slot 'FDR' must be numeric.")
    }

    TRUE
  },
  prototype = list(
    non.flat = character(),
    p.vector = numeric(0), # Empty matrix for sc.p.vector
    p.adjusted = numeric(0), # Empty numeric vector for p.adjusted
    FDR = 0 # Default FDR value is 0
  )
)
###############################################################################
#' @title Estimates
#'
#' @description
#' S4 class to store results of the polynomial term selection from full model.
#'
#' @slot sol A data frame for summary results of the stepwise regression.
#' For each selected gene the following values are given:
#'   \itemize{
#'     \item{p-value}{of the regression ANOVA.}
#'     \item{R-squared}{of the model.}
#'     \item{p-value}{of the regression coefficients of the selected variables.}
#'   }
#' @slot coefficients A data frame containing regression coefficients for the adjusted models.
#' @slot group.coeffs A matrix with the coefficients of the implicit models of each experimental group.
#' @slot t.score A data frame containing tscores for each covariate in polynomial glm.
#' @slot path A character vector containing the branching path.
#' @slot influ.info A matrix with genes containing influencial data.
#'
#' @name Estimates
#' @aliases Estimates-class
#' @rdname Estimates-class
#' @importFrom methods is new
#' @keywords classes
setClass(
  "Estimates",
  representation(
    sol = "data.frame",
    coefficients = "data.frame",
    group.coeffs = "matrix",
    t.score = "data.frame",
    path = "character",
    influ.info = "matrix"
  ),
  validity = function(object) {
    if (!is.data.frame(object@sol)) {
      stop("sol slot must be a data.frame")
    }
    if (!is.data.frame(object@coefficients)) {
      stop("coefficients slot must be a data.frame")
    }
    if (!is.matrix(object@group.coeffs)) {
      stop("group.coeffs slot must be a matrix.")
    }
    if (!is.data.frame(object@t.score)) {
      stop("t.score slot must be a data.frame")
    }
    if (!is.character(object@path)) {
      stop("path slot must be a character.")
    }
    if (!is.matrix(object@influ.info)) {
      stop("influ.info slot must be a matrix.")
    }
  },
  prototype = list(
    sol = data.frame(),
    coefficients = data.frame(),
    group.coeffs = matrix(0, nrow = 1, ncol = 1),
    t.score = data.frame(),
    path = character(),
    influ.info = matrix(NA, nrow = 0, ncol = 0)
  )
)
###############################################################################
#' @title Significant
#'
#' @description
#' S4 class to store significantly selected results features and their clusters.
#'
#' @slot sig.genes list
#' @slot feature.clusters list
#'
#' @name Significant
#' @aliases Significant-class
#' @rdname Significant-class
#' @importFrom methods is new as
#' @keywords classes
setClass(
  "Significant",
  representation(
    sig.genes = "list",
    feature.clusters = "list"
  ),
  validity = function(object) {
    if (!is.list(object@sig.genes)) {
      stop("sig.genes slot must be a list")
    }
  },
  prototype = list(
    sig.genes = list(),
    feature.clusters = list()
  )
)
###############################################################################
#' @title ScMaSigPro
#'
#' @description
#' S4 class to store analysis of ScMaSigPro worflow. It stores results,
#' parameters and associated data. Inherits slots from \code{SingleCellExperiment}.
#'
#' @slot sparse Object of Class SingleCellExperiment. See \pkg{SingleCellExperiment} for more details.
#' @slot profile Object of Class VariableProfiles See \pkg{VariableProfiles} for more details.
#' @slot estimate Object of Class Estimates See \code{\link{Estimates}} for more details.
#' @slot dense  Object of Class SingleCellExperiment. See \pkg{SingleCellExperiment} for more details.
#' @slot design Object of Class MatrixDesign. See \code{\link{MatrixDesign}} for more details.
#' @slot param Object of Class ParameterConfig. See \code{\link{ParameterConfig}} for more details.
#' @slot sig.genes ABC
#'
#' @name ScMaSigPro
#' @aliases ScMaSigPro-class
#' @rdname ScMaSigPro-class
#' @exportClass ScMaSigPro
#' @importFrom methods is new as
#' @keywords classes

setClass(
  "ScMaSigPro",
  representation(
    sparse = "SingleCellExperiment",
    profile = "VariableProfiles",
    estimate = "Estimates",
    dense = "SingleCellExperiment",
    design = "MatrixDesign",
    param = "ParameterConfig",
    sig.genes = "Significant"
  ),
  validity = function(object) {
    # Check sparse slot
    if (!validObject(object@sparse)) {
      stop("sparse slot is not a valid SingleCellExperiment object.")
    }

    # Check VariableProfiles slot
    if (!validObject(object@profile)) {
      stop("profile slot is not a valid VariableProfiles object.")
    }

    # Check Estimates slot
    if (!validObject(object@estimate)) {
      stop("Estimates slot is not a valid Estimates object.")
    }

    # Check dense slot
    if (!validObject(object@dense)) {
      stop("dense slot is not a valid SingleCellExperiment object.")
    }

    # Check MatrixDesign slot
    if (!validObject(object@design)) {
      stop("design slot is not a valid MatrixDesign object.")
    }

    # Check ParameterConfig slot
    if (!validObject(object@param)) {
      stop("param slot is not a valid ParameterConfig object.")
    }
    # Check ParameterConfig slot
    if (!validObject(object@param)) {
      stop("'sig.genes' slot is not a valid ParameterConfig object.")
    }
  },
  prototype = list(
    profile = new("VariableProfiles"), # Assuming you've defined VariableProfiles with its prototype
    estimate = new("Estimates"), # Assuming you've defined Estimates with its prototype
    param = new("ParameterConfig"), # Assuming you've defined Estimates with its prototype
    sig.genes = new("Significant")
  )
)
