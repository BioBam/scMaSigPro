"PlotGroups" <-
  function(
      data, edesign = NULL, time = edesign[, 1], groups = edesign[
        ,
        c(3:ncol(edesign))
      ], repvect = edesign[, 2], show.lines = TRUE, show.fit = FALSE,
      dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05, nvar.correction = FALSE,
      summary.mode = "median", groups.vector = NULL, main = NULL, sub = NULL,
      xlab = "Time", ylab = "Expression value", item = NULL, ylim = NULL, pch = 21,
      col = NULL, legend = TRUE, cex.legend = 1,
      lty.legend = NULL, ...) {
    if (!is.vector(data)) {
      if (summary.mode == "representative") {
        distances <- apply(as.matrix(dist(data,
          diag = TRUE,
          upper = TRUE
        )), 1, sum)
        representative <- names(distances)[distances == min(distances)]
        yy <- as.numeric(data[rownames(data) == representative, ])
        sub <- paste("Representative:", representative)
      } else if (summary.mode == "median") {
        yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
        if (is.null(sub)) {
          sub <- paste("Median profile of", nrow(data), item, sep = " ")
        }
      } else {
        stop("not valid summary.mode")
      }
      if (dim(data)[1] == 1) {
        main <- rownames(data)
        sub <- NULL
      }
    } else if (length(data) != 0) {
      yy <- as.numeric(data)
      sub <- rownames(data)
    } else {
      stop("empty data")
    }
    if (is.null(ncol(groups))) {
      ncol <- 1
      legend <- FALSE
      codeg <- "group"
    } else {
      ncol <- ncol(groups)
      codeg <- as.character(colnames(groups))
    }
    reps <- i.rank(repvect)
    y <- vector(mode = "numeric", length = length(unique(reps)))
    x <- vector(mode = "numeric", length = length(unique(reps)))
    g <- matrix(nrow = length(unique(reps)), ncol = ncol)
    for (k in 1:length(y)) {
      y[k] <- mean(yy[reps == k], na.rm = TRUE)
      x[k] <- mean(time[reps == k])
      for (j in 1:ncol) {
        g[k, j] <- mean(groups[reps == k, j])
      }
    }
    if (is.null(ylim)) {
      ylim <- c(min(as.numeric(yy), na.rm = TRUE), max(as.numeric(yy),
        na.rm = TRUE
      ))
    }
    abcissa <- x
    xlim <- c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE) *
      1.3)
    if (is.null(col)) {
      color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
      color2 <- groups
      for (j in 1:ncol) {
        color2[, j] <- color2[, j] * j
      }
      color2 <- as.vector(apply(color2, 1, sum) + 1)
    } else {
      color1 <- col
      color2 <- groups
      for (j in 1:ncol) {
        color2[, j] <- color2[, j] * col[j]
      }
      color2 <- as.vector(apply(color2, 1, sum))
    }
    plot(
      x = time, y = yy, pch = pch, xlab = xlab, ylab = ylab, main = main,
      xaxt = "n", sub = sub, ylim = ylim, xlim = xlim, col = color2, ...
    )
    arg <- list(...)
    axis(1,
      at = unique(abcissa), labels = unique(abcissa),
      lwd = 1, cex.axis = arg$cex.axis
    )
    if (show.fit) {
      rm <- matrix(yy, nrow = 1, ncol = length(yy))
      rownames(rm) <- c("ratio medio")
      colnames(rm) <- rownames(dis)
      fit.y <- T.fit(rm,
        design = dis, step.method = step.method,
        min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction
      )
      betas <- fit.y$coefficients
    }
    for (i in 1:ncol(groups)) {
      group <- g[, i]
      if ((show.fit) && !is.null(betas)) {
        li <- c(2:6)
        a <- reg.coeffs(
          coefficients = betas, groups.vector = groups.vector,
          group = colnames(groups)[i]
        )
        a <- c(a, rep(0, (7 - length(a))))
        curve(
          a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) +
            a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5),
          from = min(time), to = max(time),
          col = color1[i], add = TRUE, lty = li[i], ...
        )
      }
      if (show.lines) {
        lx <- abcissa[group != 0]
        ly <- y[group != 0]
        ord <- order(lx)
        lxo <- lx[ord]
        lyo <- ly[ord]
        lines(lxo, lyo, col = color1[i], ...)
      }
    }
    op <- par(bg = "white")
    if (legend) {
      legend(max(abcissa, na.rm = TRUE) * 1.02, ylim[1],
        legend = codeg,
        text.col = color1, col = color1, lty = lty.legend,
        yjust = 0, bty = "n", cex = cex.legend
      )
    }
    par(op)
  }
