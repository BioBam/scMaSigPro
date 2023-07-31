"average.rows" <-
  function(x, index, match, r = 0.7) {
    name <- index[match(rownames(x), as.character(match))]
    uninames <- unique(name)
    prom <- matrix(0, nrow = length(uninames), ncol = ncol(x))
    colnames(prom) <- colnames(x)
    rownames(prom) <- uninames
    for (i in 1:length(uninames)) {
      c <- x[name == uninames[i], ]
      if (is.element(TRUE, cor(t(c), use = "pairwise.complete.obs") <
        r)) {
        prom[i, ] <- c(rep(NA, ncol(c)))
        rownames(prom)[i] <- NA
      } else {
        prom[i, ] <- apply(as.matrix(c), 2, mean, na.rm = TRUE)
      }
    }
    prom <- as.data.frame(prom)
    prom <- prom[!is.na(rownames(prom)), ]
    prom
  }
