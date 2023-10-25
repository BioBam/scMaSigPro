#' Class "addParamClass"
#'
#' @slot bin_pseudotime_colname A character representing the name of the column for binned pseudotime values.
#' @slot path_prefix A character representing the prefix for path labeling.
#' @slot root_label A character representing the label for the root of the tree.
#' @slot pseudotime_colname A character representing the name of the column for pseudotime values.
#' @slot binning A character representing the binning method.
#' @slot bin_method A character representing the algorithm used for binning.
#' @slot path_colname A character representing the name of the column for path values.
#' @slot bin_colname A character representing the name of the column for bin values.
#' @slot bin_size_colname A character representing the name of the column for bin sizes.
#' @slot bin_members_colname A character representing the name of the column for bin members.
#' @slot Q Significance Level
#' @slot min.obs Minimum value to estimate the model (degree+1) x Groups + 1. (Default = 6).
#' @slot g Integer. Number of genes taken in the regression fit.
#' @slot MT.adjust Pvalue correction
#' @slot epsilon convergence tolerance
#'
#' @name addParamClass
#' @aliases addParamClass-class
#' @rdname addParamClass-class
#' @exportClass addParamClass
#' @importFrom methods is new
#' @keywords classes

setClass(
  "addParamClass",
  representation(
    bin_pseudotime_colname = "character",
    path_prefix = "character",
    root_label = "character",
    pseudotime_colname = "character",
    binning = "character",
    bin_method = "character",
    path_colname = "character",
    bin_colname = "character",
    bin_size_colname = "character",
    bin_members_colname = "character",
    g = "integer",
    Q = "numeric", # Significance level (default is 0.05)
    min.obs = "numeric", # Minimum value to estimate the model (degree+1) x Groups + 1
    MT.adjust = "character",
    epsilon = "numeric"
  ),
  validity = function(object) {
    errors <- character(0)

    # Check if any character slots are empty or not of type character
    char_slots <- c(
      "bin_pseudotime_colname", "path_prefix", "root_label",
      "pseudotime_colname", "binning", "bin_method",
      "path_colname", "bin_colname", "bin_size_colname",
      "bin_members_colname", "MT.adjust"
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
    
    # Check for slot min.obs
    if (!is.numeric(object@min.obs)) {
        stop("Slot 'min.obs' must be an integer.")
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
    bin_pseudotime_colname = "scmp_binned_pseudotime",
    path_prefix = "Path",
    root_label = "root",
    pseudotime_colname = "Pseudotime",
    binning = "universal",
    path_colname = "Path",
    bin_method = "Sturges",
    bin_colname = "scmp_bin",
    g = 0L, # Default g value is 0
    Q = 0.05,
    min.obs = 6, # Default min.obs value is 0
    bin_size_colname = "scmp_bin_size",
    bin_members_colname = "scmp_bin_members",
    MT.adjust = "BH",
    epsilon = 0.00001
  )
)
