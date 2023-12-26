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
#' for binned Pseudotime values in 'Dense' data. See `colData` from the
#' \pkg{SingleCellExperiment} package.
#' @slot path_prefix Update Description..
#' @slot root_label Update Description..
#' @slot ptime_col A character string representing the column name
#' for inferred Pseudotime values in 'Sparse' data. See `colData` from the
#' \pkg{SingleCellExperiment} package.
#' @slot bin_method A character string representing the algorithm used for
#' binning.
#' @slot path_col A character string representing the column name for branching
#' path assignment in 'Sparse' or 'Dense'data. See `colData` from the
#' \pkg{SingleCellExperiment} package.
#' @slot bin_col A character string representing the name of the column in which
#' bin labels are stored.
#' @slot bin_size_col A character string representing the name of the column in
#' which bin sizes per bin are stored.
#' @slot bin_mem_col A character string representing the name of the column in
#' which cells per bin are stored.
#' @slot anno_col A character string representing the column name for cell level
#' metadata containing cell level annotations. (Default is "cell_type").
#' @slot g Update Description..
#' @slot p_value Significance Level.
#' @slot min_na Minimum values needed per gene across cells to estimate the
#' model.
#' @slot mt_correction A character string specifying the p-value correction
#' method.
#' @slot epsilon Model convergence tolerance.
#' @slot selection_method A character string specifying the method for stepwise
#' regression.
#' @slot offset A logical value specifying whether to use offset during fitting.
#' @slot log_offset A logical value specifying whether to take the logarithm of
#' the offsets.
#' @slot max_it Maximum number of iterations to fit the model.
#' @slot poly_degree Order of the polynomial linear model.
#' @slot distribution Distribution of the error term.
#' @slot cluster_method Clustering method used for clustring significant genes.
#' @slot use_dim Dimension to use for filling the missing values before
#' clustering.
#' @slot fill_na Method to fill the missing values.
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
    ptime_col = "character",
    selection_method = "character",
    bin_method = "character",
    path_col = "character",
    bin_col = "character",
    bin_size_col = "character",
    bin_mem_col = "character",
    anno_col = "character",
    g = "integer",
    p_value = "numeric",
    min_na = "numeric",
    mt_correction = "character",
    epsilon = "numeric",
    offset = "logical",
    log_offset = "logical",
    max_it = "integer",
    poly_degree = "integer",
    distribution = "ANY",
    cluster_method = "character",
    use_dim = "character",
    fill_na = "character"
  ),
  validity = function(object) {
    errors <- character(0)

    # Check if any character slots are empty or not of type character
    char_slots <- c(
      "bin_ptime_col", "path_prefix", "root_label", "ptime_col", "bin_method",
      "path_col", "bin_col", "bin_size_col", "bin_mem_col",
      "mt_correction", "selection_method", "anno_col", "cluster_method",
      "use_dim", "fill_na"
    )

    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name)
      if (length(slot_value) == 0 || !is.character(slot_value)) {
        errors <- c(errors, paste(slot_name, "must not be empty and should
                                  be of type character."))
      }
    }

    if (!is.integer(object@g)) {
      stop("Slot 'g' must be an integer.")
    }

    if (!is.numeric(object@p_value)) {
    }

    if (!is.numeric(object@epsilon)) {
      stop("Slot 'epsilon' must be numeric.")
    }

    # Check for slot min.na
    if (!is.numeric(object@min_na)) {
      stop("Slot 'min_na' must be an integer.")
    }

    # Check if any of the character slots have multiple values
    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name)
      if (length(slot_value) > 1) {
        errors <- c(errors, paste(
          slot_name,
          "should not contain multiple values."
        ))
      }
    }

    if (length(errors) == 0) TRUE else errors
  },
  prototype = list(
    bin_ptime_col = "scmp_binned_pseudotime",
    path_prefix = "Path",
    root_label = "root",
    ptime_col = "Pseudotime",
    path_col = "Path",
    bin_method = "Sturges",
    bin_col = "scmp_bin",
    g = 0L,
    p_value = 0.05,
    min_na = 6,
    bin_size_col = "scmp_bin_size",
    bin_mem_col = "scmp_bin_members",
    mt_correction = "BH",
    selection_method = "backward",
    epsilon = 1e-8,
    offset = TRUE,
    log_offset = FALSE,
    max_it = 100L,
    anno_col = "cell_type",
    poly_degree = 2L,
    distribution = MASS::negative.binomial(10),
    cluster_method = "hclust",
    use_dim = "col",
    fill_na = "zero"
  )
)

###############################################################################
#' @title MatrixDesign
#'
#' @description
#' S4 class to store Branching Path Assignments and Prediction Matrix
#'
#' @slot predictor_matrix A matrix containing independent variables for
#' model fitting
#' @slot groups.vector A character vector specifying the branching path for each
#' term of the polynomial GLM.
#' @slot assignment_matrix A matrix containing binary assignment to branching paths
#' Additionally has two columns for assignment of binned Pseudotime and
#' replicate.
#'
#' @name MatrixDesign
#' @aliases MatrixDesign-class
#' @rdname MatrixDesign-class
#' @importFrom methods is new
#' @keywords classes
setClass(
  "MatrixDesign",
  representation(
    predictor_matrix = "matrix",
    groups.vector = "character",
    assignment_matrix = "matrix"
  ),
  prototype = list(
    predictor_matrix = matrix(NA, nrow = 0, ncol = 0),
    groups.vector = character(),
    assignment_matrix = matrix(NA, nrow = 0, ncol = 0)
  ),
  validity = function(object) {
    if (!validObject(object@predictor_matrix)) {
      stop("predictor_matrix slot is not a valid matrix.")
    }
    if (!validObject(object@groups.vector)) {
      stop("groups.vector slot is not a valid character vector.")
    }
    if (!validObject(object@assignment_matrix)) {
      stop("assignment_matrix slot is not a valid matrix.")
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
#' @slot non_flat A character vector of gene names marked as non-flat profiles.
#' @slot p_values A numeric vector containing the computed p-values for
#' non-flat profiles.
#' @slot adj_p_values A numeric vector containing the computed adjusted
#' p-values for non-flat profiles.
#' @slot fdr False Discovery Rate (FDR) value for the given significance level.
#'
#' @name VariableProfiles
#' @aliases VariableProfiles-class
#' @rdname VariableProfiles-class
#' @keywords classes

#' @author Priyansh Srivastava <spriyansh29@@gmail.com>
#' @importFrom stats family gaussian poisson
#' @importFrom utils data combn
#' @importFrom MASS negative.binomial

# Define the VariableProfiles with the following slots:
setClass("VariableProfiles",
  slots = c(
    non_flat = "character",
    p_values = "numeric",
    adj_p_values = "numeric",
    fdr = "numeric"
  ),
  validity = function(object) {
    if (!is.character(object@non_flat)) {
      stop("Slot 'non_flat' must be a character vector of gene names.")
    }
    if (!is.numeric(object@p_values)) {
      stop("Slot 'p_values' must be a numeric vector.")
    }
    if (!is.numeric(object@adj_p_values)) {
      stop("Slot 'adj_p_values' must be a numeric vector.")
    }
    if (!is.numeric(object@fdr)) {
      stop("Slot 'fdr' must be numeric.")
    }
    TRUE
  },
  prototype = list(
    non_flat = character(),
    p_values = numeric(0),
    adj_p_values = numeric(0),
    fdr = numeric(0)
  )
)
###############################################################################
#' @title Estimates
#'
#' @description
#' S4 class to store results of the polynomial term selection from full model.
#'
#' @slot significance_matrix A matrix storing p-values for each gene per
#' polynomial term, R2 and global p-values.
#' @slot coefficient_matrix A matrix storing beta estimates for each polynomial
#' term.
#' @slot path_coefficient_matrix A matrix with the coefficients of each
#' branching path.
#' @slot t_score_matrix A matrix storing the t-scores for each polynomial term.
#' @slot path A character vector containing the branching path.
#' @slot influential A matrix with genes containing influential data.
#'
#' @name Estimates
#' @aliases Estimates-class
#' @rdname Estimates-class
#' @importFrom methods is new
#' @keywords classes
setClass(
  "Estimates",
  representation(
    significance_matrix = "matrix",
    coefficient_matrix = "matrix",
    path_coefficient_matrix = "matrix",
    t_score_matrix = "matrix",
    path = "character",
    influential = "matrix"
  ),
  validity = function(object) {
    if (!is.matrix(object@significance_matrix)) {
      stop("significance_matrix slot must be a matrix.")
    }
    if (!is.matrix(object@coefficient_matrix)) {
      stop("coefficient_matrix slot must be a matrix.")
    }
    if (!is.matrix(object@path_coefficient_matrix)) {
      stop("path_coefficient_matrix slot must be a matrix.")
    }
    if (!is.matrix(object@t_score_matrix)) {
      stop("t_score_matrix slot must be a matrix.")
    }
    if (!is.character(object@path)) {
      stop("path slot must be a character.")
    }
    if (!is.matrix(object@influential)) {
      stop("influential slot must be a matrix.")
    }
  },
  prototype = list(
    significance_matrix = matrix(NA, nrow = 0, ncol = 0),
    coefficient_matrix = matrix(NA, nrow = 0, ncol = 0),
    path_coefficient_matrix = matrix(0, nrow = 0, ncol = 0),
    t_score_matrix = matrix(NA, nrow = 0, ncol = 0),
    path = character(),
    influential = matrix(NA, nrow = 0, ncol = 0)
  )
)
###############################################################################
#' @title Significant
#'
#' @description
#' S4 class to store significantly selected features and their clusters.
#'
#' @slot genes Named list of genes selected based on user-specified thresholds.
#' @slot clusters Named list of clusters selected based on user-specified
#' thresholds.
#'
#' @name Significant
#' @aliases Significant-class
#' @rdname Significant-class
#' @importFrom methods is new as
#' @keywords classes
setClass(
  "Significant",
  representation(
    genes = "list",
    clusters = "list"
  ),
  validity = function(object) {
    if (!is.list(object@genes)) {
      stop("genes slot must be a list")
    }
    if (!is.list(object@clusters)) {
      stop("clusters slot must be a list")
    }
  },
  prototype = list(
    genes = list(),
    clusters = list()
  )
)
###############################################################################
#' @title ScMaSigPro
#'
#' @description
#' S4 class to store analysis of ScMaSigPro worflow. It stores results,
#' parameters and associated data. Inherits slots from \code{SingleCellExperiment}.
#'
#' @slot Sparse Object of class \code{\link{SingleCellExperiment}}.
#' See \pkg{SingleCellExperiment} for more details.
#' @slot Dense Object of class \code{\link{SingleCellExperiment}}.
#' See \pkg{SingleCellExperiment} for more details.
#' @slot Design Object of Class \code{\link{MatrixDesign}}.
#' @slot Profile Object of Class \code{\link{VariableProfiles}}.
#' @slot Estimate Object of Class \code{\link{Estimates}}.
#' @slot Significant Object of Class \code{\link{Significant}}.
#' @slot Parameters Object of Class \code{\link{ParameterConfig}}.
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
    Sparse = "SingleCellExperiment",
    Dense = "SingleCellExperiment",
    Design = "MatrixDesign",
    Profile = "VariableProfiles",
    Estimate = "Estimates",
    Significant = "Significant",
    Parameters = "ParameterConfig"
  ),
  validity = function(object) {
    # Check sparse slot
    if (!validObject(object@Sparse)) {
      stop("'Sparse' slot is not a valid 'SingleCellExperiment' class object.")
    }

    # Check VariableProfiles slot
    if (!validObject(object@Profile)) {
      stop("'Profile' slot is not a valid 'VariableProfiles' class object.")
    }

    # Check Estimates slot
    if (!validObject(object@Estimate)) {
      stop("'Estimates' slot is not a valid 'Estimates' class object.")
    }

    # Check dense slot
    if (!validObject(object@Dense)) {
      stop("'Dense' slot is not a valid 'SingleCellExperiment' class object.")
    }

    # Check MatrixDesign slot
    if (!validObject(object@Design)) {
      stop("'Design' slot is not a valid 'MatrixDesign' class object.")
    }

    # Check ParameterConfig slot
    if (!validObject(object@Parameters)) {
      stop("'Parameters' slot is not a valid 'ParameterConfig' class object.")
    }
    # Check ParameterConfig slot
    if (!validObject(object@Significant)) {
      stop("'Significant' slot is not a valid 'Significant' object.")
    }
  },
  prototype = list(
    Profile = new("VariableProfiles"),
    Estimate = new("Estimates"),
    Parameters = new("ParameterConfig"),
    Significant = new("Significant")
  )
)
