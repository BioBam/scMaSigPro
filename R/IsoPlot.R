IsoPlot <- function(get, name, only.sig.iso = FALSE, ylim = NULL, xlab = "Time", ylab = "Expression value", points = TRUE, cex.main = 3, cex.legend = 1.5) {
  edesign <- get$Model$design$edesign
  data <- get$Model$data
  gen <- get$Model$gen

  if (only.sig.iso) {
    sig.iso2 <- get$get2$summary
    gen <- as.character(gen[rownames(data) %in% sig.iso2])
    data <- get$get2$sig.genes$sig.profiles
  }

  data <- data[gen == name, ]
  nt <- nrow(data)
  time <- edesign[, 1]
  groups <- edesign[, c(3:ncol(edesign))]
  repvect <- edesign[, 2]
  n.groups <- ncol(edesign) - 2

  if (n.groups == 1) {
    group <- factor(groups, labels = colnames(edesign)[3])
  } else {
    group <- factor(as.matrix(edesign[, 3:ncol(edesign)]) %*% c(1:n.groups), labels = colnames(edesign)[3:ncol(edesign)])
  }

  yy <- as.numeric(data[1, ])

  if (n.groups == 1) {
    ncol <- 1
    codeg <- colnames(edesign)[3]
  } else {
    ncol <- ncol(groups)
    codeg <- as.character(colnames(groups))
  }
  if (is.null(ylim)) ylims <- c(min(data, na.rm = TRUE), max(data, na.rm = TRUE))

  #----------------------------------------------------------------------
  # POINTS-PLOT FOR EACH GROUP:
  #----------------------------------------------------------------------

  op <- c(5, 4, 4, 0)
  yaxt <- NULL
  par(mfrow = c(1, ncol + 1))

  for (i in 1:ncol)
  {
    if (ncol == 2 && i == 2) {
      yaxt <- "n"
      op <- c(5, 0, 4, 4)
    }

    if (ncol > 2) {
      if (i == ncol) {
        yaxt <- "n"
        op <- c(5, 0, 4, 4)
      }
      if (i > 1 && i < ncol) {
        yaxt <- "n"
        op <- c(5, 2, 4, 2)
      }
    }
    ploti(dati = data[, group == codeg[i]], ti = time[group == codeg[i]], op = op, ylim = ylims, ylab = ylab, xlab = xlab, yaxt = yaxt, main = codeg[i], points, cex.main)
  }
  plot.new()
  legend("left", legend = rownames(data), title = " ", cex = cex.legend, lty = 1, col = c(2:(nt + 1)), bty = "n", lwd = 2)
  legend("left", legend = rep(" ", nt), title = name, cex = cex.main, bty = "n", title.adj = 1)
}


#-------------------------------------------------------------------------------------
# AUXILIAR FUNCTION AUXILIAR ploti()
#-------------------------------------------------------------------------------------
ploti <- function(dati, ti, op, ylim = ylim, ylab = ylab, xlab = xlab, yaxt = yaxt, main = main, points, cex.main) {
  nt <- nrow(dati)
  y1 <- as.numeric(dati[1, ])
  par(mar = op)

  tii <- tapply(ti, ti, mean)
  plot(tii, tapply(y1, ti, mean), type = "l", ylim = ylim, col = 2, xaxt = "n", ylab = ylab, xlab = xlab, yaxt = yaxt, main = main, cex.main = cex.main)
  axis(1, at = unique(ti), labels = unique(ti))

  if (points) {
    points(ti, y1, col = 2)
  }

  for (j in 2:nt) {
    yj <- as.numeric(dati[j, ])
    lines(tii, tapply(yj, ti, mean), col = (j + 1))
    if (points) {
      points(ti, yj, ylim = ylim, col = (j + 1))
    }
  }
}
