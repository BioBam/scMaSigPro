"PlotProfiles" <-
  function(
      data, cond, cex.axis = 0.5, ylim = NULL,
      repvect, main = NULL, sub = NULL, color.mode = "rainbow", item = NULL) {
    pos.vline <- repvect[1:length(repvect)] - c(0, repvect[1:(length(repvect) -
      1)])
    for (i in 1:length(pos.vline)) {
      if (pos.vline[i] != 0) {
        pos.vline[i] <- i
      } else {
        pos.vline[i] <- pos.vline[i - 1]
      }
    }
    if (is.null(ylim)) {
      ylim <- c(min(as.matrix(data), na.rm = TRUE) * 1.1, max(as.matrix(data),
        na.rm = TRUE
      ) * 1.1)
    }
    if (!is.vector(data)) {
      n <- dim(data)[2]
      m <- dim(data)[1]
      if (m == 1) {
        nom <- rownames(data)
      } else {
        nom <- paste("Cluster", main, "(", m, item, ")", sep = " ")
      }
      plot(
        x = c(1:n), y = data[1, ], type = "l", col = 1,
        ylim = ylim, ylab = "expression value", xlab = " ",
        main = nom, xaxt = "n"
      )
      axis(1,
        at = 1:n, labels = substr(cond, 1, 26), cex.axis = cex.axis,
        las = 2
      )
      if (color.mode == "rainbow") {
        abline(v = pos.vline, col = "light gray")
        for (i in 1:dim(data)[1]) {
          lines(x = c(1:n), y = data[i, ], col = i)
        }
      } else if (color.mode == "gray") {
        for (i in 1:dim(data)[1]) {
          lines(x = c(1:n), y = data[i, ], col = "light gray")
        }
        yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
        lines(x = c(1:n), y = yy, col = "black")
      } else {
        stop("Invalid mode, must be one of rainbow, gray")
      }
    } else {
      n <- length(data)
      plot(
        x = c(1:n), y = data, type = "l", col = 1, ylim = ylim,
        ylab = "expression value", sub, xaxt = "n", xlab = " "
      )
      axis(1,
        at = 1:n, labels = cond, cex.axis = cex.axis,
        las = 2
      )
      abline(v = pos.vline, col = "light gray")
    }
  }
