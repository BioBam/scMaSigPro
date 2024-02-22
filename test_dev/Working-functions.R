# functions

# log(0)
log.cond <- function(x) {
  sol <- numeric(length(x))
  for (i in 1:length(x)) {
    if (x[i] == 0) {
      sol[i] <- 0
    } else {
      sol[i] <- log(x[i])
    }
  }
  sol
}

# R2 from deviance in glm
R2 <- function(model) {
  (model$null.deviance - model$deviance) / model$null.deviance
}

#---------------------------------------------------------------------------
# plot-function to use several times:
plot.pred <- function(y, model, main, type = "response", ylab = "y") {
  plot(Bin, y, col = as.numeric(Bpath) + 1, main = main, ylab = ylab)
  pred <- predict(model, type = type)
  lines(Bin[1:8], pred[1:8], col = 2)
  lines(Bin[9:16], pred[9:16], col = 3)
}
