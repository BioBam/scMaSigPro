get.siggenes <- function(
    tstep, rsq = 0.7, add.IDs = FALSE, IDs = NULL, matchID.col = 1,
    only.names = FALSE, vars = c("all", "each", "groups"), significant.intercept = "dummy",
    groups.vector = NULL, trat.repl.spots = "none", index = IDs[
      ,
      (matchID.col + 1)
    ], match = IDs[, matchID.col], r = 0.7) {
  dis <- tstep$dis
  edesign <- tstep$edesign
  groups.vector <- tstep$groups.vector
  not.all.empty <- function(x) (is.element(FALSE, x == " "))
  indep <- strsplit(
    colnames(dis)[grep("x", colnames(dis))[1]],
    "x"
  )[[1]][1]
  if (any(tstep$sol[, 2] > rsq)) {
    sig.pvalues <- tstep$sol[which(tstep$sol[, 2] > rsq), ]
    sig.profiles <- tstep$sig.profiles[which(tstep$sol[
      ,
      2
    ] > rsq), ]
    coefficients <- tstep$coefficients[which(tstep$sol[
      ,
      2
    ] > rsq), ]
    group.coeffs <- tstep$group.coeffs[which(tstep$sol[
      ,
      2
    ] > rsq), ]
    if (vars == "all") {
      sigs <- sig.profiles
      summary <- rownames(sig.profiles)
      if (only.names) {
        sigs <- rownames(sig.profiles)
      }
      if (add.IDs) {
        row.names <- rownames(sig.profiles)
        ids <- IDs[match(rownames(sig.profiles), IDs[
          ,
          matchID.col
        ]), ]
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
        sigs, coeffs, gc, ps, nrow(sig.profiles),
        edesign, groups.vector
      )
      names(sig.genes) <- c(
        "sig.profiles", "coefficients",
        "group.coeffs", "sig.pvalues", "g", "edesign",
        "groups.vector"
      )
    } else if (vars == "each") {
      sig.genes <- as.list(paste("var", c("independ", colnames(dis)),
        sep = "."
      ))
      summary <- matrix(" ", ncol = ncol(dis) + 1, nrow = nrow(sig.profiles))
      colnames(summary) <- c("independ", colnames(dis))
      for (i in 1:ncol(summary)) {
        sigs <- sig.profiles[which(!is.na(sig.pvalues[
          ,
          (2 + i)
        ])), ]
        coeffs <- coefficients[which(!is.na(sig.pvalues[
          ,
          (2 + i)
        ])), ]
        gc <- group.coeffs[which(!is.na(sig.pvalues[
          ,
          (2 + i)
        ])), ]
        ps <- sig.pvalues[which(!is.na(sig.pvalues[
          ,
          (2 + i)
        ])), ]
        if (nrow(sigs) > 0) {
          names.sigs <- rownames(sigs)
        } else {
          names.sigs <- NULL
        }
        summary[, i] <- c(names.sigs, rep(" ", nrow(sig.profiles) -
          nrow(sigs)))
        sig.genes[[i]] <- list(
          sigs, coeffs, gc, ps,
          nrow(sigs), edesign, groups.vector
        )
        names(sig.genes[[i]]) <- c(
          "sig.profiles", "coefficients",
          "group.coeffs", "sig.pvalues", "g", "edesign",
          "groups.vector"
        )
      }
      names(sig.genes) <- c("independ", colnames(dis))
      summary <- as.data.frame(summary[apply(
        summary, 1,
        not.all.empty
      ), ])
    } else if (vars == "groups") {
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
        # In case 1 series, degree=1 and significant intercept="dummy", only time variable is in group.sig:
        if (length(selc) == 1) {
          colnames(group.sig) <- cnames1[selc]
        }

        cnames2 <- colnames(group.sig)
        group.sig <- as.data.frame(group.sig[, groups.vector[selc] ==
          group[i]])
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
          sigs, coeffs, gc, ps,
          nrow(ps), edesign, groups.vector
        )
        names(sig.genes[[i]]) <- c(
          "sig.profiles", "coefficients",
          "group.coeffs", "sig.pvalues", "g", "edesign",
          "groups.vector"
        )
      }
      names(sig.genes) <- unique(groups.vector)
      if (nrow(summary) > 1) {
        summary <- as.data.frame(summary[apply(
          summary,
          1, not.all.empty
        ), ])
      }
    } else {
      stop("invalid vars value, must be one of: all, each, groups")
    }
    if (trat.repl.spots == "average") {
      if (vars != "all") {
        for (i in 1:length(sig.genes)) {
          sig.genes[[i]][[1]] <- average.rows(sig.genes[[i]][[1]],
            index = index, match = match, r = r
          )
          for (j in c(2:3)) {
            sig.genes[[i]][[j]] <- average.rows(sig.genes[[i]][[j]],
              index = index, match = match, r = -1
            )
            sig.genes[[i]][[j]] <- sig.genes[[i]][[j]][is.element(
              sig.genes[[i]][[j]],
              sig.genes[[i]][[1]]
            ), ]
          }
        }
      } else {
        sig.genes[[1]] <- average.rows(sig.genes[[1]],
          index = index, match = match, r = r
        )
        for (j in c(2:3)) {
          sig.genes[[j]] <- average.rows(sig.genes[[j]],
            index = index, match = match, r = -1
          )
          sig.genes[[j]] <- sig.genes[[j]][is.element(
            sig.genes[[j]],
            sig.genes[[1]]
          ), ]
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
            ids <- IDs[match(
              rownames(sig.genes2[[i]][[1]]),
              IDs[, matchID.col]
            ), ]
          } else {
            stop("function parameters no compatible (add.IDs, trat.repl.spots)")
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
  output <- list(sig.genes, summary)
  names(output) <- c("sig.genes", "summary")
  output
}
