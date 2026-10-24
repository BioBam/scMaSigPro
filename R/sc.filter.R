#' @title List the significant genes
#'
#' @description
#' `sc.filter()` lists the genes whose selected model ([sc.t.fit()]) has an
#' R-squared of at least `rsq` and a p-value of at most `p_value`: in one list,
#' one list per term or one list per path (`vars`), as `get.siggenes()` of
#' maSigPro does.
#'
#' @param rsq The smallest R-squared of the selected model of a gene.
#' @param p_value The largest p-value of the selected model of a gene; by
#'   default that of [sc.p.vector()].
#' @param term_p_value The largest p-value of a term, for `vars = "each"` and
#'   `vars = "groups"`.
#' @param vars `"all"`: one list with every gene that passes; `"each"`: one
#'   list per term, with the genes whose term passes `term_p_value`;
#'   `"groups"`: one list per path, with the genes with a term of that path
#'   that passes it.
#' @param intercept With `vars = "groups"`, the intercept terms that count:
#'   `"dummy"`, the dummies of the paths (the differences of their intercepts
#'   from the reference); `"all"`, these and the intercept of the reference
#'   path; `"none"`, neither (with more than two paths, only the dummy of the
#'   first path after the reference is left out).
#'
#' @return The object with the gene lists in the `Significant` slot
#' (class [Significant-class]); earlier clusters are removed.
#'
#' @seealso [plotIntersect()] for the overlaps of the lists;
#' [Significant-class].
#'
#' @references
#' Conesa A, Nueda MJ, Ferrer A, Talón M (2006). maSigPro: a method to
#' identify significantly differential expression profiles in time-course
#' microarray experiments. Bioinformatics 22(9), 1096-1102.
#' \doi{10.1093/bioinformatics/btl056}
#'
#' Nueda MJ, Tarazona S, Conesa A (2014). Next maSigPro: updating maSigPro
#' bioconductor package for RNA-seq time series. Bioinformatics 30(18),
#' 2598-2602. \doi{10.1093/bioinformatics/btu333}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#'
#' @family workflow steps
#' @inheritParams scmp-params
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' scmp <- sc.filter(scmp.ob, rsq = 0.6, vars = "groups")
#' res <- scmpResults(scmp)
#' # The number of genes in each list
#' c(Path1 = sum(res$sig_Path1), Path2vsPath1 = sum(res$sig_Path2vsPath1))
#'
#' @export
sc.filter <- function(scmpObj,
                      rsq = 0.7,
                      p_value = scmpObj@Parameters@p_value,
                      vars = c("all", "each", "groups"),
                      intercept = "dummy",
                      term_p_value = 0.05,
                      includeInflu = TRUE) {
    # Check Validity of the object
    .check_scmp(scmpObj)

    # The default lists the choices and means the first, as in match.arg()
    if (missing(vars)) {
        vars <- "all"
    }
    vars <- .check_choice(vars, c("all", "each", "groups"),
        msg = "Invalid selection for 'vars'"
    )
    if (vars == "groups") {
        intercept <- .check_choice(intercept, c("dummy", "none", "all"),
            msg = "Invalid selection for 'intercept'"
        )
    }

    # Initiate empty list
    sig.list <- list()

    # Extract sol
    sol <- showSol(scmpObj, includeInflu = includeInflu)

    # Why: a gene without a p-value or R-squared must fail the
    # thresholds below.
    # Replace the NAs with corresponding values
    sol[is.na(sol$`p-value`), "p-value"] <- 1
    sol[is.na(sol$`R-squared`), "R-squared"] <- 0
    sol[, !(colnames(sol) %in% c("p-value", "R-squared"))][is.na(sol[, !(colnames(sol) %in% c("p-value", "R-squared"))])] <- 1

    # Filter based on RSQ and P-value
    sol <- sol[sol$`p-value` <= p_value, ]
    sol <- sol[sol$`R-squared` >= rsq, ]

    # If`all`
    if (vars == "all") {
        # Return a list of gene names
        sig.list[["all"]] <- rownames(sol)
    } else if (vars == "each") {
        # Subset
        sol.sub <- sol[, !(colnames(sol) %in% c("p-value", "R-squared")), drop = FALSE]

        # Get the genes of each term: lapply(), because apply() returns a
        # matrix when every term selects the same number of genes
        sig.list <- lapply(sol.sub, function(column) {
            row.names(sol.sub)[which(column <= term_p_value)]
        })
        names(sig.list) <- sub("p.valor_", "", names(sig.list), fixed = TRUE)
    } else if (vars == "groups") {
        # Subset
        sol.sub <- sol[, !(colnames(sol) %in% c("p-value", "R-squared")), drop = FALSE]

        # Get group_vector, from T fit
        group_vector <- scmpObj@Estimate@path

        # Based on the dummy, none and all
        if (intercept == "all") {
            select_cols <- seq(from = 1, to = ncol(sol.sub))
        } else if (intercept == "dummy") {
            select_cols <- seq(from = 2, to = ncol(sol.sub))
        } else if (intercept == "none") {
            select_cols <- seq(from = 3, to = ncol(sol.sub))
        }

        # subset groupVector and sol
        group_vector <- group_vector[select_cols]
        sol.sub <- sol.sub[, select_cols, drop = FALSE]

        # Create Named list for traversal
        group_term_list <- as.list(as.vector(colnames(sol.sub)))
        names(group_term_list) <- group_vector

        # Available groups
        avail_groups <- unique(group_vector)

        # Get per path
        sig.list <- lapply(avail_groups, function(group_i) {
            # Get names per path
            group_term_list_sub <- group_term_list[names(group_term_list) %in% group_i]

            # get subset
            sol.sub.term <- sol.sub[, unlist(group_term_list_sub), drop = FALSE]

            # Get group wise genes
            sig.list.tmp <- row.names(sol.sub.term)[apply(sol.sub.term, 1, function(x) any(x <= term_p_value))]

            # Return
            return(sig.list.tmp)
        })

        # Add names
        names(sig.list) <- avail_groups
    }

    # Create Object
    siggenes.object <- new("Significant",
        genes = sig.list
    )

    # Update the slot
    scmpObj@Significant <- siggenes.object
    validObject(scmpObj)

    # Return
    return(scmpObj)
}
