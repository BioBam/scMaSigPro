# ScMaSigProClass with all the super-classes
# Author: Priyansh Srivastava

###############################################################################
#' @title ParameterConfig
#'
#' @description
#' The settings of an analysis: the names of the columns of the cell and bin
#' metadata and the parameters of every step of the workflow.
#'
#' @slot bin_ptime_col The column of the bin metadata with the binned
#'   pseudotime.
#' @slot path_prefix The default path column of [as_scmp()].
#' @slot root_label The label of the root cells given to [as_scmp()], kept for
#'   reference.
#' @slot ptime_col The column of the cell metadata with the pseudotime.
#' @slot bin_method The rule for the number of bins of [sc.squeeze()].
#' @slot path_col The column of the cell and bin metadata with the path
#'   labels.
#' @slot bin_col The column of the bin metadata with the bin names.
#' @slot bin_size_col The column of the bin metadata with the number of
#'   cells of each bin.
#' @slot bin_mem_col The column of the bin metadata with the cells of each
#'   bin.
#' @slot anno_col The column of the cell metadata with cell annotations; not
#'   used by the steps.
#' @slot g The number of genes that [sc.p.vector()] fitted: those left after
#'   removing the genes with fewer than `min_na` values, with zero counts in
#'   every bin, or below `min_counts` or `min_bins_expressed`.
#' @slot p_value The level of [sc.p.vector()] and [sc.t.fit()].
#' @slot min_na The minimum number of bins with a value (not `NA`) that
#'   [sc.p.vector()] required; it only removes genes from counts with missing
#'   values.
#' @slot mt_correction The multiple-testing correction of [sc.p.vector()].
#' @slot link The link function of the fits.
#' @slot epsilon The convergence tolerance of the fits.
#' @slot selection_method The stepwise method of [sc.t.fit()].
#' @slot offset Whether the fits used the offset of each bin, the logarithm of
#'   its size factor (see [sc.p.vector()]).
#' @slot max_it The maximum number of iterations of the fits.
#' @slot poly_degree The degree of the polynomial of [sc.set.poly()].
#' @slot distribution The family of the fits, a [stats::family()] object;
#'   with an estimated theta it also holds `theta_estimate` and `theta_gene`
#'   (see [sc.p.vector()]).
#' @slot cluster_method The clustering method of [sc.cluster.trend()].
#' @slot use_dim How [sc.cluster.trend()] filled missing values: by column or
#'   by row.
#' @slot fill_na What [sc.cluster.trend()] filled missing values with.
#' @slot min_counts The minimum total count over the bins that
#'   [sc.p.vector()] required of a gene (0: no minimum).
#' @slot min_bins_expressed The minimum number of bins with a count above 0
#'   that [sc.p.vector()] required of a gene (0: no minimum).
#' @slot test The test of the full model of each gene that [sc.p.vector()]
#'   was asked for: `"chisq"`, the chi-squared test, or `"F"`, the F test.
#'
#' @details
#' `validObject()` checks that every slot except `distribution` holds a
#' single non-missing value, that `min_counts` and `min_bins_expressed` are
#' at least 0, that `test` is `"chisq"` or `"F"` and that `distribution` is a
#' [stats::family()] object.
#'
#' @name ParameterConfig
#' @aliases ParameterConfig-class
#' @return An object of class `ParameterConfig`: the settings of an
#' analysis, one value per slot. Every step of the workflow records the
#' settings it used in the `Parameters` slot of a [ScMaSigPro-class] object,
#' and [showParams()] lists them.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' showParams(scmp.ob)
#'
#' @rdname ParameterConfig-class
#' @family classes
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
        max_it = "integer",
        poly_degree = "integer",
        distribution = "ANY",
        cluster_method = "character",
        use_dim = "character",
        fill_na = "character",
        min_counts = "numeric",
        min_bins_expressed = "numeric",
        test = "character"
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
        for (s in c("min_counts", "min_bins_expressed")) {
            value <- slot(object, s)
            if (length(value) == 1L && !is.na(value) && value < 0) {
                msg <- c(msg, paste0("'", s, "' must be at least 0"))
            }
        }
        test <- object@test
        known <- test %in% c("chisq", "F")
        if (length(test) == 1L && !is.na(test) && !known) {
            msg <- c(msg, "'test' must be \"chisq\" or \"F\"")
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
        max_it = 100L,
        anno_col = "cell_type",
        poly_degree = 2L,
        distribution = MASS::negative.binomial(10),
        cluster_method = "hclust",
        use_dim = "col",
        fill_na = "zero",
        min_counts = 0,
        min_bins_expressed = 0,
        test = "chisq"
    )
)

###############################################################################
#' @title MatrixDesign
#'
#' @description
#' The design of the polynomial model of [sc.set.poly()].
#'
#' @slot predictor_matrix The design matrix: one row per bin, one column per
#'   term ([predictors()]).
#' @slot groups.vector The path of each column of `predictor_matrix`: the
#'   reference path for its powers of pseudotime, `<path>vs<ref>` for the
#'   terms of another path.
#' @slot assignment_matrix The path assignment matrix: one row per bin, with
#'   the binned pseudotime, a replicate index and a 0/1 column per path
#'   ([pathAssign()]).
#' @slot offset The offset of each bin, the logarithm of its size factor (see
#'   [sc.p.vector()]), or an empty vector when the fits used no offset or no
#'   fit was made yet.
#'
#' @details
#' `validObject()` checks that `assignment_matrix` has the rows of
#' `predictor_matrix` (one per bin), that `offset` is empty or has one value
#' per row and that `groups.vector` has one value per column of
#' `predictor_matrix`.
#'
#' @name MatrixDesign
#' @aliases MatrixDesign-class
#' @return An object of class `MatrixDesign`: the design of the polynomial
#' model, which [sc.set.poly()] stores in the `Design` slot of a
#' [ScMaSigPro-class] object, with the offsets that [sc.p.vector()] adds.
#' [predictors()] and [pathAssign()] return its matrices.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' head(predictors(scmp.ob))
#'
#' @rdname MatrixDesign-class
#' @family classes
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
        n_offset <- length(object@offset)
        if (n_offset != 0L && n_offset != nrow(pred)) {
            msg <- c(msg, paste(
                "'offset' must be empty or have one value per row of",
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
#' The test of the full model of each gene by [sc.p.vector()].
#'
#' @slot non_flat The genes with a non-flat profile.
#' @slot p_values The p-value of every fitted gene, named by gene.
#' @slot adj_p_values The adjusted p-values, in the same order.
#' @slot fdr The largest p-value of the genes with a non-flat profile, empty
#'   when there is none.
#'
#' @details
#' `validObject()` checks that `p_values` and `adj_p_values` have the same
#' length and names, that `fdr` holds at most one value and that every gene
#' in `non_flat` is named in `p_values`.
#'
#' @name VariableProfiles
#' @aliases VariableProfiles-class
#' @return An object of class `VariableProfiles`: the test of the full model
#' of each gene, which [sc.p.vector()] stores in the `Profile` slot of a
#' [ScMaSigPro-class] object. [scmpResults()] returns its p-values.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' res <- scmpResults(scmp.ob)
#' head(res[, c("p_value", "adj_p_value", "non_flat")])
#'
#' @rdname VariableProfiles-class
#' @family classes
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
#' The models that [sc.t.fit()] selected.
#'
#' @slot significance_matrix For each gene with a significant selected
#'   model, its p-value and R-squared and the p-value of each term
#'   ([showSol()]).
#' @slot coefficient_matrix The coefficients ([showCoeff()]).
#' @slot path_coefficient_matrix The coefficients of the polynomial of each
#'   path ([showGroupCoeff()]).
#' @slot t_score_matrix The t statistics ([showTS()]).
#' @slot path The path of the intercept and of each term.
#' @slot influential The genes with influential observations
#'   ([showInflu()]).
#'
#' @details
#' `validObject()` checks that `t_score_matrix` has the shape and rows of
#' `coefficient_matrix`, that `path_coefficient_matrix` has its rows, and
#' that every row of `significance_matrix` is a row of `coefficient_matrix`.
#'
#' @name Estimates
#' @aliases Estimates-class
#' @return An object of class `Estimates`: the selected model of each gene,
#' which [sc.t.fit()] stores in the `Estimate` slot of a [ScMaSigPro-class]
#' object. [showSol()] and the other result readers return its matrices.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' head(showCoeff(scmp.ob))
#'
#' @rdname Estimates-class
#' @family classes
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
#' The gene lists of [sc.filter()] and the clusters of [sc.cluster.trend()].
#'
#' @slot genes A named list of gene lists.
#' @slot clusters A list with the cluster of each clustered gene, named by
#'   gene.
#'
#' @details
#' `validObject()` checks that every element of `genes` is a character
#' vector of gene names and that `clusters` holds one value per gene, named
#' by gene.
#'
#' @name Significant
#' @aliases Significant-class
#' @return An object of class `Significant`: the gene lists of
#' [sc.filter()] and the clusters of [sc.cluster.trend()], stored in the
#' `Significant` slot of a [ScMaSigPro-class] object. [scmpResults()]
#' returns them per gene.
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' table(scmpResults(scmp.ob)$cluster)
#'
#' @rdname Significant-class
#' @family classes
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
# The version of scMaSigPro that is running, as a string such as "0.0.5".
# getNamespaceVersion() is base R (utils is not imported) and reads the
# loaded namespace, so it is the version of the code that runs.
.scmp_version <- function() unname(getNamespaceVersion("scMaSigPro"))

###############################################################################
#' @title ScMaSigPro
#'
#' @description
#' An analysis of the scMaSigPro workflow: the data at the level of the cells
#' and of the pseudotime bins, as two `SingleCellExperiment` objects, and the
#' results and settings of every step.
#'
#' @slot Sparse The cells: a `SingleCellExperiment` with their counts and
#'   metadata, including the pseudotime and path columns.
#' @slot Dense The pseudotime bins: a `SingleCellExperiment` with their
#'   pseudo-bulk counts and metadata.
#' @slot Design The design, a [MatrixDesign-class] object.
#' @slot Profile The test of each gene, a [VariableProfiles-class] object.
#' @slot Estimate The selected models, an [Estimates-class] object.
#' @slot Significant The gene lists and clusters, a [Significant-class]
#'   object.
#' @slot Parameters The settings, a [ParameterConfig-class] object.
#' @slot version The version of scMaSigPro that created the object or last
#'   changed its structure in
#'   \code{\link[=updateObject,ScMaSigPro-method]{updateObject()}}.
#'
#' @details
#' `validObject()` checks every slot with the validity method of its class
#' and reports each problem with the name of the slot; `version` must be a
#' single version number. The constructors, the steps of the workflow,
#' [pb_counts()] and the setters check the object they return, so a step or
#' setter that would return an invalid object stops with these problems
#' instead. Objects saved by an earlier version of the package, for example
#' with `saveRDS()`, are brought up to date with
#' \code{\link[=updateObject,ScMaSigPro-method]{updateObject()}}.
#'
#' @section Why a class of its own:
#' `ScMaSigPro` holds two `SingleCellExperiment` objects because the analysis
#' works at two resolutions whose columns differ. `Sparse` holds the cells,
#' with their pseudotime and path labels. `Dense` holds the pseudotime bins
#' (the pseudo-bulk) to which the polynomial GLMs are fitted. One
#' `SingleCellExperiment` cannot hold both: `altExp()` needs the columns of
#' the main experiment, and bins stored in `metadata()` would lose the
#' validity checks and accessors of an experiment. The other slots hold
#' results that have no place in the SummarizedExperiment model: the design
#' matrices (`MatrixDesign`), the per-gene fits of the full model
#' (`VariableProfiles`), the coefficients and significance of the selected
#' models (`Estimates`), the gene lists and clusters (`Significant`) and a
#' record of every parameter used (`ParameterConfig`), for reproducibility.
#' Both experiments are built and changed through the SummarizedExperiment
#' accessors only, and `validObject()` runs after every step of the
#' workflow. To use the results with other tools, [scmpResults()] returns
#' them as a `DataFrame` with one row per gene, and [scmpToRowData()] writes
#' them into the `rowData()` of a `SingleCellExperiment`.
#'
#' @name ScMaSigPro
#' @aliases ScMaSigPro-class
#' @return An object of class `ScMaSigPro`. [as_scmp()] and [create_scmp()]
#' make one with the cells in the `Sparse` slot; [sc.squeeze()] adds the
#' bins (`Dense`), [sc.set.poly()] the design (`Design`), [sc.p.vector()]
#' the test of each gene (`Profile`), [sc.t.fit()] the selected models
#' (`Estimate`), [sc.filter()] and [sc.cluster.trend()] the gene lists and
#' clusters (`Significant`), and every step its settings (`Parameters`).
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' scmp.ob
#' validObject(scmp.ob)
#'
#' @rdname ScMaSigPro-class
#' @family classes
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
        Parameters = "ParameterConfig",
        version = "character"
    ),
    validity = function(object) {
        msg <- character()
        # By default validObject() checks only the class of each slot.
        for (s in setdiff(slotNames(object), "version")) {
            problems <- validObject(slot(object, s), test = TRUE)
            if (!isTRUE(problems)) {
                msg <- c(msg, paste0("slot '", s, "': ", problems))
            }
        }
        v <- object@version
        one_version <- length(v) == 1L && !is.na(v) &&
            !is.na(package_version(v, strict = FALSE))
        if (!one_version) {
            msg <- c(msg, "'version' must be a single version number")
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
        Significant = new("Significant"),
        # Set when the package is installed or loaded, so every new object
        # carries the version of the code that creates it.
        version = .scmp_version()
    )
)
