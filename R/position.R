"position" <-
  function(matrix, vari) {
    a <- colnames(matrix)
    b <- a == vari
    c <- c(1:length(a))
    d <- c[b]
    return(d)
  }
