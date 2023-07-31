"suma2Venn" <- function(x, size = 30, cexil = 0.9, cexsn = 1, zcolor = heat.colors(ncol(x)), ...) {
  G <- ncol(x)
  L <- vector("list", G)
  names(L) <- colnames(x)
  for (i in 1:G)
  {
    y <- as.character(x[, i])
    y <- y[y != " "]
    L[[i]] <- y
  }

  venn(L, size = size, cexil = cexil, cexsn = cexsn, zcolor = zcolor, ...)
}
