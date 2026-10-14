# ScMaSigProClass with all the super-classes
# Author: Priyansh Srivastava

###############################################################################
#' @title ParameterConfig
#'
#' @description
#' S4 Class to hold parameters used in the analysis. It is constantly updated
#' as the analysis progresses.
#'
#' @slot bin_ptime_col A character string representing the column name
#' for binned Pseudotime values in 'Dense' data. See `colData` from the
#' \pkg{SingleCellExperiment} package.
#' @slot path_prefix A character string, the prefix for branching path
#' labels given to [as_scmp()] (Default is "Path"); [as_scmp()] also uses
#' it as the default `path_col`.
#' @slot root_label A character string, the label for root cells given to
#' [as_scmp()] (Default is "root"). It is stored for reference only.
#' @slot ptime_col A character string representing the column name
#' for inferred Pseudotime values in 'Sparse' data. See `colData` from the
#' \pkg{SingleCellExperiment} package.
#' @slot bin_method A character string representing the algorithm used for
#' binning.
#' @slot path_col A character string representing the column name for branching
#' path assignment in 'Sparse' or 'Dense' data. See `colData` from the
#' \pkg{SingleCellExperiment} package.
#' @slot bin_col A character string representing the name of the column in which
#' bin labels are stored.
#' @slot bin_size_col A character string representing the name of the column in
#' which bin sizes per bin are stored.
#' @slot bin_mem_col A character string representing the name of the column in
#' which cells per bin are stored.
#' @slot anno_col A character string representing the column name for cell level
#' metadata containing cell level annotations. (Default is "cell_type").
#' @slot g An integer, the number of features that [sc.p.vector()] fitted:
#' those left after removing features with fewer than `min_na` values or with
#' zero counts in every bin.
#' @slot p_value Significance Level.
#' @slot min_na Minimum values needed per gene across cells to estimate the
#' model.
#' @slot mt_correction A character string specifying the p-value correction
#' method.
#' @slot link Type of link function to use in the model. Default is "log".
#' @slot epsilon Model convergence tolerance.
#' @slot selection_method A character string specifying the method for stepwise
#' regression.
#' @slot offset A logical value specifying whether to use offset during fitting.
#' @slot log_offset A logical value specifying whether to take the logarithm of
#' the offsets.
#' @slot max_it Maximum number of iterations to fit the model.
#' @slot poly_degree Order of the polynomial linear model.
#' @slot distribution Distribution of the error term.
#' @slot cluster_method Clustering method used for clustering significant genes.
#' @slot use_dim Dimension to use for filling the missing values before
#' clustering.
#' @slot fill_na Method to fill the missing values.
#'
#' @details
#' `validObject()` checks that every slot except `distribution` holds a
#' single non-missing value and that `distribution` is a [stats::family()]
#' object.
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
        link = "character",
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
        msg <- character()
        # S4 already enforces each slot's type; the workflow also needs one
        # value per slot.
        for (s in setdiff(slotNames(object), "distribution")) {
            value <- slot(object, s)
            if (length(value) != 1L || is.na(value)) {
                msg <- c(msg, paste0(
                    "'", s, "' must be a single non-missing value"
                ))
            }
        }
        if (!inherits(object@distribution, "family")) {
            msg <- c(msg, "'distribution' must be a 'family' object")
        }
        if (length(msg)) msg else TRUE
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
        link = "log",
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
#' @slot offset A numeric vector containing the offset values.
#'
#' @details
#' `validObject()` checks that `assignment_matrix` has the rows of
#' `predictor_matrix` (one per bin), that `offset` has one value per row and
#' that `groups.vector` has one value per column of `predictor_matrix`.
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
        assignment_matrix = "matrix",
        offset = "numeric"
    ),
    prototype = list(
        predictor_matrix = matrix(NA, nrow = 0, ncol = 0),
        groups.vector = character(),
        assignment_matrix = matrix(NA, nrow = 0, ncol = 0),
        offset = numeric()
    ),
    validity = function(object) {
        msg <- character()
        pred <- object@predictor_matrix
        alloc <- object@assignment_matrix
        same_rows <- nrow(alloc) == nrow(pred) &&
            identical(rownames(alloc), rownames(pred))
        if (!same_rows) {
            msg <- c(msg, paste(
                "'assignment_matrix' must have the rows of",
                "'predictor_matrix'"
            ))
        }
        if (length(object@offset) != nrow(pred)) {
            msg <- c(msg, paste(
                "'offset' must have one value per row of",
                "'predictor_matrix'"
            ))
        }
        if (length(object@groups.vector) != ncol(pred)) {
            msg <- c(msg, paste(
                "'groups.vector' must have one value per column of",
                "'predictor_matrix'"
            ))
        }
        if (length(msg)) msg else TRUE
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
#' @details
#' `validObject()` checks that `p_values` and `adj_p_values` have the same
#' length and names, that `fdr` holds at most one value and that every gene
#' in `non_flat` is named in `p_values`.
#'
#' @name VariableProfiles
#' @aliases VariableProfiles-class
#' @rdname VariableProfiles-class
#' @keywords classes

#' @author Priyansh Srivastava <spriyansh29@@gmail.com>
#' @importFrom stats gaussian
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
        msg <- character()
        p <- object@p_values
        adj <- object@adj_p_values
        if (length(adj) != length(p) || !identical(names(adj), names(p))) {
            msg <- c(msg, paste(
                "'adj_p_values' must have the length and names of",
                "'p_values'"
            ))
        }
        if (length(object@fdr) > 1L) {
            msg <- c(msg, "'fdr' must hold at most one value")
        }
        if (!all(object@non_flat %in% names(p))) {
            msg <- c(msg, paste(
                "every gene in 'non_flat' must be named in",
                "'p_values'"
            ))
        }
        if (length(msg)) msg else TRUE
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
#' @details
#' `validObject()` checks that `t_score_matrix` has the shape and rows of
#' `coefficient_matrix`, that `path_coefficient_matrix` has its rows, and
#' that every row of `significance_matrix` is a row of `coefficient_matrix`.
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
        msg <- character()
        genes <- rownames(object@coefficient_matrix)
        tscore <- object@t_score_matrix
        same_shape <- identical(dim(tscore), dim(object@coefficient_matrix)) &&
            identical(rownames(tscore), genes)
        if (!same_shape) {
            msg <- c(msg, paste(
                "'t_score_matrix' must have the shape and rows of",
                "'coefficient_matrix'"
            ))
        }
        if (!identical(rownames(object@path_coefficient_matrix), genes)) {
            msg <- c(msg, paste(
                "'path_coefficient_matrix' must have the rows of",
                "'coefficient_matrix'"
            ))
        }
        if (!all(rownames(object@significance_matrix) %in% genes)) {
            msg <- c(msg, paste(
                "rows of 'significance_matrix' must be rows of",
                "'coefficient_matrix'"
            ))
        }
        if (length(msg)) msg else TRUE
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
#' @details
#' `validObject()` checks that every element of `genes` is a character
#' vector of gene names and that `clusters` holds one value per gene, named
#' by gene.
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
        msg <- character()
        if (!all(vapply(object@genes, is.character, logical(1)))) {
            msg <- c(msg, paste(
                "every element of 'genes' must be a character vector of",
                "gene names"
            ))
        }
        cl <- object@clusters
        named <- !is.null(names(cl)) && all(nzchar(names(cl)))
        if (length(cl) > 0L && (!named || any(lengths(cl) != 1L))) {
            msg <- c(msg, paste(
                "'clusters' must hold one value per gene, named",
                "by gene"
            ))
        }
        if (length(msg)) msg else TRUE
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
#' @slot Sparse Object of class \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' See \pkg{SingleCellExperiment} for more details.
#' @slot Dense Object of class \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' See \pkg{SingleCellExperiment} for more details.
#' @slot Design Object of Class \code{\link{MatrixDesign}}.
#' @slot Profile Object of Class \code{\link{VariableProfiles}}.
#' @slot Estimate Object of Class \code{\link{Estimates}}.
#' @slot Significant Object of Class \code{\link{Significant}}.
#' @slot Parameters Object of Class \code{\link{ParameterConfig}}.
#'
#' @details
#' `validObject()` checks every slot with the validity method of its class
#' and reports each problem with the name of the slot. Objects are created
#' with [create_scmp()] or [as_scmp()].
#'
#' @name ScMaSigPro
#' @aliases ScMaSigPro-class
#' @rdname ScMaSigPro-class
#' @exportClass ScMaSigPro
#' @importFrom methods is new as validObject
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
        msg <- character()
        # By default validObject() checks only the class of each slot.
        for (s in slotNames(object)) {
            problems <- validObject(slot(object, s), test = TRUE)
            if (!isTRUE(problems)) {
                msg <- c(msg, paste0("slot '", s, "': ", problems))
            }
        }
        if (length(msg)) msg else TRUE
    },
    prototype = list(
        # new("SingleCellExperiment") is not a valid object; the constructor is.
        Sparse = SingleCellExperiment(),
        Dense = SingleCellExperiment(),
        Profile = new("VariableProfiles"),
        Estimate = new("Estimates"),
        Parameters = new("ParameterConfig"),
        Significant = new("Significant")
    )
)
