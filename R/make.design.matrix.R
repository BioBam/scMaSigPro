"make.design.matrix" <-
  function(edesign, degree = 2, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign))) {
    control.label <- colnames(edesign)[group.cols][1]
    if (dim(as.matrix(edesign))[2] > 3) {
      dummy.cols <- group.cols[2:length(group.cols)]
      treatm.label <- paste(colnames(edesign)[dummy.cols],
        "vs", control.label,
        sep = ""
      )
      groups.label <- c(control.label, treatm.label)
      matrix.dummy <- as.matrix(edesign[, dummy.cols])
      ## Shared origin
      dummy <- NULL
      j <- 0
      origen <- min(edesign[, time.col])
      origen <- edesign[edesign[, 1] == origen, ]
      for (i in 1:length(dummy.cols)) {
        share <- apply(origen[, c(3, dummy.cols[i])], 1, sum)
        if (!is.element(TRUE, share > 1)) {
          j <- j + 1
          dummy <- cbind(dummy, matrix.dummy[, i])
          colnames(dummy)[j] <- treatm.label[i]
        }
      }
      time <- as.matrix(edesign[, time.col])
      colnames(time) <- colnames(edesign)[time.col]
      dis <- cbind(dummy, time)
      rownames(dis) <- rownames(edesign)
      groups.vector <- c(colnames(dummy), control.label)
      colnames.dis <- colnames(dis)
      dis <- cbind(dis, dis[, ncol(dis)] * matrix.dummy)
      colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col],
        "x", colnames(edesign)[dummy.cols],
        sep = ""
      ))
      groups.vector <- c(groups.vector, treatm.label)
      if (degree >= 2) {
        for (i in 2:degree) {
          colnames.dis <- colnames(dis)
          dis <- cbind(dis, edesign[, time.col]^i, edesign[
            ,
            time.col
          ]^i * edesign[, dummy.cols])
          colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col],
            i,
            sep = ""
          ), paste(colnames(edesign)[time.col],
            "", i, "x", colnames(edesign)[dummy.cols],
            sep = ""
          ))
          groups.vector <- c(groups.vector, groups.label)
        }
      }
    } else {
      dis <- as.matrix(edesign[, time.col])
      colnames(dis) <- colnames(edesign)[time.col]
      rownames(dis) <- rownames(edesign)
      if (degree > 1) {
        for (i in 2:degree) {
          colnames.dis <- colnames(dis)
          dis <- cbind(dis, edesign[, time.col]^i)
          colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col],
            i,
            sep = ""
          ))
        }
      }
      groups.vector <- rep(colnames(edesign)[group.cols], degree)
    }
    output <- list(dis, groups.vector, edesign)
    names(output) <- c("dis", "groups.vector", "edesign")
    output
  }
