#' Extract significant genes for sets of variables in time series gene expression experiments
#'
#' This function creates lists of significant genes for a set of variables whose significance
#' value has been computed with the \code{T.fit} function.
#'
#' @title Extract significant genes for sets of variables in time series gene expression experiments
#'
#' @description
#' There are 3 possible values for the \code{vars} argument:
#' \itemize{
#'   \item \code{"all"}: generates one single matrix or gene list with all significant genes.
#'   \item \code{"each"}: generates as many significant genes extractions as variables in the general regression model.
#'   \item \code{"groups"}: generates a significant genes extraction for each experimental group.
#' }
#'
#' @param tstep A \code{T.fit} object.
#' @param rsq Cut-off level at the R-squared value for the stepwise regression fit. Only genes with R-squared more than 'rsq' are selected.
#' @param add.IDs Logical indicating whether to include additional gene id's in the result.
#' @param IDs Matrix containing additional gene id information (required when \code{add.IDs = TRUE}).
#' @param matchID.col Number of the matching column in the matrix \code{IDs} for adding gene ids.
#' @param only.names Logical. If \code{TRUE}, expression values are omitted in the results.
#' @param vars Variables for which to extract significant genes.
#' @param significant.intercept Experimental groups for which significant intercept coefficients are considered.
#' @param groups.vector Required when \code{vars = "groups"}.
#' @param trat.repl.spots Treatment given to replicate spots. Possible values are \code{"none"} and \code{"average"}.
#' @param index Argument of the \code{\link{average.rows}} function to use when \code{trat.repl.spots = "average"}.
#' @param match Argument of the \code{\link{average.rows}} function to use when \code{trat.repl.spots = "average"}.
#' @param r Minimum Pearson correlation coefficient for replicated spots profiles to be averaged.
#'
#' @details
#' Refer to the function description for details on the arguments and their usage.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{summary}: A vector or matrix listing significant genes for the variables given by the function parameters.
#'   \item \code{sig.genes}: A list with detailed information on the significant genes found for the variables given by the function parameters.
#'   Each element of the list is also a list containing:
#'     \describe{
#'       \item{\code{sig.profiles}:}{Expression values of significant genes.}
#'       \item{\code{coefficients}:}{Regression coefficients of the adjusted models.}
#'       \item{\code{group.coeffs}:}{Regression coefficients of the implicit models of each experimental group.}
#'       \item{\code{sig.pvalues}:}{P-values of the regression coefficients for significant genes.}
#'       \item{\code{g}:}{Number of genes.}
#'       \item{\code{...}:}{Arguments passed by previous functions.}
#'     }
#'   }
#'
#' @references
#' Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. (2006).
#' maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
#' Bioinformatics, 22(9), 1096-1102. \url{https://doi.org/10.1093/bioinformatics/btl056}
#'
#' @author Ana Conesa and Maria Jose Nueda (mj.nueda@ua.es)
#'
#' @examples
#' # Example usage of the function can be placed here.
#'
#' @keywords manip
#'
#' @export
get.siggenes <- function(tstep, rsq = 0.7, add.IDs = FALSE, IDs = NULL, matchID.col = 1,
                         only.names = FALSE, vars = c("all", "each", "groups"),
                         significant.intercept = "dummy",
                         groups.vector = NULL, trat.repl.spots = "none",
                         index = IDs[, (matchID.col + 1)], match = IDs[, matchID.col],
                         r = 0.7) {
  # Extract data from the tstep object
  dis <- tstep$dis
  edesign <- tstep$edesign
  groups.vector <- tstep$groups.vector

  # Function to check if all elements of a vector are not empty
  not.all.empty <- function(x) (is.element(FALSE, x == " "))

  # Extract the independent variable from the column names of the dis data
  indep <- strsplit(colnames(dis)[grep("x", colnames(dis))[1]], "x")[[1]][1]

  # Check if there are significant genes based on the rsq cutoff
  if (any(tstep$sol[, 2] > rsq)) {
    # Filter significant genes and their associated information
    sig.pvalues <- tstep$sol[which(tstep$sol[, 2] > rsq), ]
    sig.profiles <- tstep$sig.profiles[which(tstep$sol[, 2] > rsq), ]
    coefficients <- tstep$coefficients[which(tstep$sol[, 2] > rsq), ]
    group.coeffs <- tstep$group.coeffs[which(tstep$sol[, 2] > rsq), ]

    # Case 1: Extract all significant genes together
    if (vars == "all") {
      sigs <- sig.profiles
      summary <- rownames(sig.profiles)
      if (only.names) {
        sigs <- rownames(sig.profiles)
      }
      if (add.IDs) {
        row.names <- rownames(sig.profiles)
        ids <- IDs[match(rownames(sig.profiles), IDs[, matchID.col]), ]
        if (!only.names) {
          sigs <- cbind(ids, sigs)
        } else {
          sigs <- ids
        }
        rownames(sigs) <- row.names
      }
      coeffs <- coefficients
      gc <- group.coeffs
      ps <- sig.pvalues
      sig.genes <- list(
        sigs, coeffs, gc, ps, nrow(sig.profiles), edesign, groups.vector
      )
      names(sig.genes) <- c(
        "sig.profiles", "coefficients", "group.coeffs", "sig.pvalues",
        "g", "edesign", "groups.vector"
      )
    }
    # Case 2: Extract significant genes separately for each variable
    else if (vars == "each") {
      sig.genes <- as.list(paste("var", c("independ", colnames(dis)), sep = "."))
      summary <- matrix(" ", ncol = ncol(dis) + 1, nrow = nrow(sig.profiles))
      colnames(summary) <- c("independ", colnames(dis))
      for (i in 1:ncol(summary)) {
        sigs <- sig.profiles[which(!is.na(sig.pvalues[, (2 + i)])), ]
        coeffs <- coefficients[which(!is.na(sig.pvalues[, (2 + i)])), ]
        gc <- group.coeffs[which(!is.na(sig.pvalues[, (2 + i)])), ]
        ps <- sig.pvalues[which(!is.na(sig.pvalues[, (2 + i)])), ]
        if (nrow(sigs) > 0) {
          names.sigs <- rownames(sigs)
        } else {
          names.sigs <- NULL
        }
        summary[, i] <- c(names.sigs, rep(" ", nrow(sig.profiles) - nrow(sigs)))
        sig.genes[[i]] <- list(
          sigs, coeffs, gc, ps, nrow(sigs), edesign, groups.vector
        )
        names(sig.genes[[i]]) <- c(
          "sig.profiles", "coefficients", "group.coeffs", "sig.pvalues", "g",
          "edesign", "groups.vector"
        )
      }
      names(sig.genes) <- c("independ", colnames(dis))
      summary <- as.data.frame(summary[apply(summary, 1, not.all.empty), ])
    }
    # Case 3: Extract significant genes for each experimental group
    else if (vars == "groups") {
      if (is.null(groups.vector)) {
        if (is.null(tstep$groups.vector)) {
          stop("groups.vector is missing")
        } else {
          groups.vector <- tstep$groups.vector
        }
      }
      group <- unique(groups.vector)
      summary <- matrix(" ", ncol = length(group), nrow = nrow(sig.profiles))
      colnames(summary) <- group
      sig.genes <- as.list(group)

      # Select coefficients based on the significant intercept argument
      if (significant.intercept == "all") {
        selc <- c(1:length(groups.vector))
      } else if (significant.intercept == "dummy") {
        selc <- c(2:length(groups.vector))
      } else if (significant.intercept == "none") {
        selc <- grep(indep, colnames(tstep$coefficients))
      } else {
        stop("invalid significant.intercept value, must be one of: all, dummy, none")
      }

      for (i in 1:length(group)) {
        group.sig <- sig.pvalues[, grep("p.valor", colnames(sig.pvalues))]
        cnames1 <- colnames(group.sig)
        group.sig <- as.data.frame(group.sig[, selc])

        # In case of 1 series, degree=1, and significant intercept="dummy", only time variable is in group.sig
        if (length(selc) == 1) {
          colnames(group.sig) <- cnames1[selc]
        }

        cnames2 <- colnames(group.sig)
        group.sig <- as.data.frame(group.sig[, groups.vector[selc] == group[i]])
        if (ncol(group.sig) == 1) {
          group.sig <- as.data.frame(group.sig)
          colnames(group.sig) <- cnames2[groups.vector[selc] == group[i]]
        }
        ps <- sig.pvalues[which(apply(group.sig, 1, function(x) {
          any(!is.na(x))
        })), ]
        sigs <- sig.profiles[which(apply(group.sig, 1, function(x) {
          any(!is.na(x))
        })), ]
        coeffs <- coefficients[which(apply(group.sig, 1, function(x) {
          any(!is.na(x))
        })), ]
        gc <- group.coeffs[which(apply(group.sig, 1, function(x) {
          any(!is.na(x))
        })), ]
        if (nrow(sigs) > 0) {
          names.sigs <- rownames(sigs)
        } else {
          names.sigs <- NULL
        }
        summary[, i] <- c(names.sigs, rep(" ", nrow(sig.profiles) - nrow(sigs)))
        sig.genes[[i]] <- list(
          sigs, coeffs, gc, ps, nrow(ps), edesign, groups.vector
        )
        names(sig.genes[[i]]) <- c(
          "sig.profiles", "coefficients", "group.coeffs", "sig.pvalues", "g",
          "edesign", "groups.vector"
        )
      }
      names(sig.genes) <- unique(groups.vector)
      if (nrow(summary) > 1) {
        summary <- as.data.frame(summary[apply(summary, 1, not.all.empty), ])
      }
    } else {
      stop("invalid vars value, must be one of: all, each, groups")
    }

    # Check if treatment replicate spots are to be averaged
    if (trat.repl.spots == "average") {
      if (vars != "all") {
        for (i in 1:length(sig.genes)) {
          sig.genes[[i]][[1]] <- average.rows(sig.genes[[i]][[1]], index = index, match = match, r = r)
          for (j in c(2:3)) {
            sig.genes[[i]][[j]] <- average.rows(sig.genes[[i]][[j]], index = index, match = match, r = -1)
            sig.genes[[i]][[j]] <- sig.genes[[i]][[j]][is.element(sig.genes[[i]][[j]], sig.genes[[i]][[1]]), ]
          }
        }
      } else {
        sig.genes[[1]] <- average.rows(sig.genes[[1]], index = index, match = match, r = r)
        for (j in c(2:3)) {
          sig.genes[[j]] <- average.rows(sig.genes[[j]], index = index, match = match, r = -1)
          sig.genes[[j]] <- sig.genes[[j]][is.element(sig.genes[[j]], sig.genes[[1]]), ]
        }
      }
    }

    sig.genes2 <- sig.genes
    if (only.names && vars != "all") {
      for (i in 1:length(sig.genes)) {
        if (!is.null(dim(sig.genes[[i]][[1]]))) {
          sig.genes[[i]][[1]] <- rownames(sig.genes[[i]][[1]])
        }
      }
    }
    if (add.IDs && vars != "all") {
      for (i in 1:length(sig.genes)) {
        if (nrow(sig.genes2[[i]][[1]]) > 1) {
          row.names <- rownames(sig.genes2[[i]][[1]])
          if (trat.repl.spots == "none") {
            ids <- IDs[match(rownames(sig.genes2[[i]][[1]]), IDs[, matchID.col]), ]
          } else {
            stop("function parameters not compatible (add.IDs, trat.repl.spots)")
          }
          if (!only.names) {
            sig.genes[[i]][[1]] <- cbind(ids, sig.genes[[i]][[1]])
          } else {
            sig.genes[[i]][[1]] <- ids
          }
          rownames(sig.genes[[i]][[1]]) <- row.names
        }
      }
    }
  } else {
    sig.genes <- NULL
    summary <- c("no significant genes")
    print("no significant genes")
  }

  # Create the output list
  output <- list(sig.genes, summary)
  names(output) <- c("sig.genes", "summary")
  output
}
