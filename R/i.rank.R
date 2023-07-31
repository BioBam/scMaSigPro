"i.rank" <-
  function(x) {
    xx <- x
    for (i in 1:length(xx)) {
      xx[i] <- c(1:length(unique(x)))[x[i] == unique(x)]
    }
    xx
  }
