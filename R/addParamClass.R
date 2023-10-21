#' @title Class "addParamClass"
#'
#' @description A class for additional parameters related to pseudotime and path analyses.
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
#'
#' @name scmpAddParamClass
#' @aliases addParam-class
#' @rdname addParamClass-class
#' @exportClass addParamClass
#' @keywords classes
#' @export
#' 
#' @importFrom methods is new
#'
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
    bin_members_colname = "character"
  ),
  validity = function(object) {
    errors <- character(0)

    # Check if any character slots are empty or not of type character
    char_slots <- c(
      "bin_pseudotime_colname", "path_prefix", "root_label",
      "pseudotime_colname", "binning", "bin_method",
      "path_colname", "bin_colname", "bin_size_colname",
      "bin_members_colname"
    )

    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name) # Corrected line
      if (length(slot_value) == 0 || !is.character(slot_value)) { # Corrected line
        errors <- c(errors, paste(slot_name, "must not be empty and should be of type character."))
      }
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
    bin_size_colname = "scmp_bin_size",
    bin_members_colname = "scmp_bin_members"
  )
)
