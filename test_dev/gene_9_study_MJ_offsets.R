## Load Packages

source("test_dev/Working-functions.R")
library(scMaSigPro)
library(MASS)
library(tidyverse)

#----------------------------------------------------------------------
## study of Gen9
#----------------------------------------------------------------------
data("scmp.ob", package = "scMaSigPro")
dim(eSparse(scmp.ob)) # eSparse es para data raw
dim(eDense(scmp.ob)) # eDense es para acumulados
data <- eSparse(scmp.ob)
Bdata <- eDense(scmp.ob)

gen9 <- data[9, ]
Bgen9 <- Bdata[9, ]

names(cSparse(scmp.ob))
Pseudotime <- cSparse(scmp.ob)$Pseudotime
names(Pseudotime) <- rownames(cSparse(scmp.ob))
BinnedPseudotime <- cSparse(scmp.ob)$scmp_binned_pseudotime
names(BinnedPseudotime) <- rownames(cSparse(scmp.ob))
Path <- cSparse(scmp.ob)$Path

# datos agrupados:
Bin <- cDense(scmp.ob)$scmp_binned_pseudotime
names(Bin) <- rownames(cDense(scmp.ob))
Bpath <- rep(c(1:2), each = 8)

table(BinnedPseudotime, Path) # cells in each bin-path
# check sums:
sum(gen9) == sum(Bgen9)
sum1 <- sum2 <- NULL
for (i in 1:8)
{
  sum1 <- c(sum1, sum(gen9[BinnedPseudotime == i]))
  sum2 <- c(sum2, sum(Bgen9[Bin == i]))
}
sum1
sum2
#----------------------------------------------------------------------
## QUESTION 1: Why sum1 is not equal to sum2?
#----------------------------------------------------------------------

# I checked in detail what is happening. The issue lies within the `eSparse()`.
# This function is a generic and it internally uses the `SummarizedExperiment:::assay()`.
# When we extract the counts using `eSparse()` it meses up the names of the cells.
# Maybe there is an update in the `SummarizedExperiment` package and they must have
# changed something with the indexing. I will check this in detail, for now I would
# depricate the eSparse() and eDense().

# I show the problem here

# Extracting the data using the eSparse and eDense functions
data <- eSparse(scmp.ob)
Bdata <- eDense(scmp.ob)

data_gen9 <- data[9, ]
Bdata_Bgen9 <- Bdata[9, ]

# Extracting the data directly from the slot
slot_gen9 <- as.matrix(scmp.ob@Sparse@assays@data@listData$counts)[9, ]
slot_Bgen9 <- as.matrix(scmp.ob@Dense@assays@data@listData$bulk.counts)[9, ]

# Check Whether the data is the same
all(data_gen9 == slot_gen9)
all(Bgen9 == slot_Bgen9)

# Check Whether the name is same
all(names(data_gen9) == names(slot_gen9)) ## This is the problem
all(names(Bgen9) == names(slot_Bgen9))

# Now let's check using the slot functions
Pseudotime <- cSparse(scmp.ob)$Pseudotime
names(Pseudotime) <- rownames(cSparse(scmp.ob))
BinnedPseudotime <- cSparse(scmp.ob)$scmp_binned_pseudotime
names(BinnedPseudotime) <- rownames(cSparse(scmp.ob))
Path <- cSparse(scmp.ob)$Path

# datos agrupados:
Bin <- cDense(scmp.ob)$scmp_binned_pseudotime
names(Bin) <- rownames(cDense(scmp.ob))
Bpath <- rep(c(1:2), each = 8)

table(BinnedPseudotime, Path) # cells in each bin-path
# check sums:
sum(slot_gen9) == sum(slot_Bgen9)
sum1_new <- sum2_new <- NULL
for (i in 1:8) {
  cell_vector <- names(BinnedPseudotime[BinnedPseudotime == i])
  cell_count_vector <- slot_gen9[names(slot_gen9) %in% cell_vector]
  sum1_new <- c(sum1_new, sum(cell_count_vector))

  bin_vector <- names(Bin[Bin == i])
  bin_count_vector <- slot_Bgen9[names(Bin) %in% bin_vector]
  sum2_new <- c(sum2_new, sum(bin_count_vector))
}
all(sum1_new == sum2_new)
# This is TRUE

# This is not an issue as within the package, as I always acess data directly from the slots.
# The eSparse and eDense are just for the user to access the data easily.

#######################################################################
# study offsets
#######################################################################
MGeo <- function(x) {
  exp(mean(log.cond(x)))
}
mGenes <- apply(Bdata, 1, MGeo)
apply(Bdata / mGenes, 2, median)

# This is not your offset_values. I suppose for the glm we must apply -1:
computedOffset <- apply(Bdata / mGenes, 2, median) - 1
offset_values <- scmp.ob@Design@offset
computedOffset - offset_values
#----------------------------------------------------------------------
## QUESTION 2: Why they are not exactly the same?
#----------------------------------------------------------------------

# The difference is due to the fact that the offset values are calculated using the
# log(estimateSizeFactorsForMatrix(counts)) function. This function is approach is
# in the DESeq2 package. Here is the function
scmp_estimateSizeFactorsForMatrix <- function(counts, locfunc = stats::median,
                                              geoMeans, controlGenes,
                                              type = c("ratio", "poscounts")) {
  type <- match.arg(type, c("ratio", "poscounts"))
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    if (type == "ratio") {
      loggeomeans <- MatrixGenerics::rowMeans(log(counts))
    } else if (type == "poscounts") {
      lc <- log(counts)
      lc[!is.finite(lc)] <- 0
      loggeomeans <- MatrixGenerics::rowMeans(lc)
      allZero <- MatrixGenerics::rowSums(counts) == 0
      loggeomeans[allZero] <- -Inf
    }
  } else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if (!(is.numeric(controlGenes) | is.logical(controlGenes))) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes, , drop = FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  if (incomingGeoMeans) {
    # stabilize size factors to have geometric mean of 1
    sf <- sf / exp(mean(log(sf)))
  }
  sf
}

# Also the offset values you are comparing against they are not the same, because
# when the we run the tsetp function they are updated. The input then is basically
# the dense matrix with the significant non-flat profiles only.

## Here is how to verify it

# Get data
slot_Bdata <- scmp.ob@Dense@assays@data@listData$bulk.counts

## Since the object is S4 we can easily update the slots
scmp.ob <- sc.p.vector(scmp.ob, offset = T, log_offset = T)

## Get the offset values
offset_values_after_pvector <- as.numeric(scmp.ob@Design@offset)

## Calculate the offset manually
offset_values_after_pvector_manual <- as.numeric(log(scmp_estimateSizeFactorsForMatrix(
  counts = slot_Bdata + 1
)))

# Check
offset_values_after_pvector - offset_values_after_pvector_manual
# This is ~0

# I use testthat framework which is used CRAN nowadays
testthat::expect_equal(
  expected = offset_values_after_pvector,
  object = offset_values_after_pvector_manual
)

######################################################################
# 4 plots that summarize the process
par(mfrow = c(2, 2))
plot(Pseudotime, gen9, col = as.numeric(Path) + 1)
title("Gen9 - row data")
plot(BinnedPseudotime, gen9, col = as.numeric(Path) + 1)
title("Gen9 - row data")
plot(Bin, Bgen9, col = Bpath + 1)
lines(Bin[1:8], Bgen9[1:8], col = 2)
lines(Bin[9:16], Bgen9[9:16], col = 3)
title("Gen9 - binned-data")

plot(Bin, log(Bgen9), col = Bpath + 1)
lines(Bin[1:8], log(Bgen9[1:8]), col = 2)
lines(Bin[9:16], log(Bgen9[9:16]), col = 3)
title("Gen9 - log(binned-data)")

######################################################################
# Model without offsets (as I did before)
#---------------------------------------------------------------------
mod0 <- glm(Bgen9 ~ Bin * as.factor(Bpath), family = negative.binomial(10))
R2(mod0)
plot.pred(y = Bgen9, model = mod0, main = "Raw data,no offset. R2=0.2706")


#--------------------------------------------------------------
# including offsets in the model you are normalizing data through the model.
# offset values are different for each data-point:
plot(offset_values, log(Bgen9))
cor(offset_values, Bgen9) # high correlation

# Comparing raw data with this a fitted model with offsets doesn't show reality
# high R2 indicates data is near the fitted values. It doesn't depends on the scale of the plot
# because is a relative measure

#------------------------------------------------------------------------------
# Model comparison with offsets (only first step, it is enough to show this)
#------------------------------------------------------------------------------
mod1 <- glm(Bgen9 ~ Bin * as.factor(Bpath), family = negative.binomial(10), offset = offset_values)
summary(mod1)
R2(mod1)

y <- Bgen9 / exp(offset_values)
mod2 <- glm(y ~ Bin * as.factor(Bpath), family = negative.binomial(10))
summary(mod2)
R2(mod2) # they are not exactly the same, but very similar

par(mfrow = c(2, 2))
plot.pred(y = Bgen9, model = mod0, main = "Raw data,no offset. R2=0.2706")
plot.pred(y = Bgen9, model = mod1, main = "Raw data-model with offset.R2=0.9174")
plot.pred(y = y, model = mod2, main = "Corrected data.R2=0.9174")
plot.pred(y = log(y), model = mod2, type = "link", main = "Corrected log(data).R2=0.9174", ylab = "log(y)")

# if you show high R2 with first graph is confusing.
plotTrend(scmp.ob, "Gene9", logs = T, logType = "log")
#----------------------------------------------------------------------
## QUESTION 3: Don't you think graphs with corrected data will be more informative?
#----------------------------------------------------------------------
