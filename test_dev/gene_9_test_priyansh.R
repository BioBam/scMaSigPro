## study of Gen9

## Load Packages
library(scMaSigPro)
library(MASS)

data("scmp.ob", package = "scMaSigPro")
dim(eSparse(scmp.ob))
dim(eDense(scmp.ob))

data <- eSparse(scmp.ob)
Bdata <- eDense(scmp.ob)
gen9 <- data[9, ]
Bgen9 <- Bdata[9, ]

names(cSparse(scmp.ob))
Pseudotime <- cSparse(scmp.ob)$Pseudotime
names(Pseudotime) <- rownames(cSparse(scmp.ob))
BinnedPseudotime <- cDense(scmp.ob)$scmp_binned_pseudotime
names(BinnedPseudotime) <- rownames(cDense(scmp.ob))
Path <- cSparse(scmp.ob)$Path

par(mfrow = c(2, 2))
plot(Pseudotime, gen9, col = as.numeric(Path) + 1)
title("Gen9 - row data")
plot(BinnedPseudotime, Bgen9, col = as.numeric(Path) + 1)
title("Gen9 - row data")

Bin <- rep(c(1:8), 2)
Bpath <- rep(c(1:2), each = 8)
plot(Bin, Bgen9, col = Bpath + 1)
lines(Bin[1:8], Bgen9[1:8], col = 2)
lines(Bin[9:16], Bgen9[9:16], col = 3)
title("Gen9 - binned-data")

plot(Bin, log(Bgen9), col = Bpath + 1)
lines(Bin[1:8], log(Bgen9[1:8]), col = 2)
lines(Bin[9:16], log(Bgen9[9:16]), col = 3)
title("Gen9 - log(binned-data). R2=0.2706")

mod <- glm(Bgen9 ~ Bin * as.factor(Bpath), family = negative.binomial(10))
summary(mod)
mod0 <- glm(Bgen9 ~ 1, family = negative.binomial(10))
anova(mod0, mod, test = "Chisq")


R2 <- function(model) {
  (model$null.deviance - model$deviance) / model$null.deviance
}
R2(mod)

#######################################################################
# Computing the R2 with offset

# Offsets are Stored in the object
# They are for the samples
offset_values <- scmp.ob@Design@offset
offset_values

# Model with offset
mod_w_offset <- glm(Bgen9 ~ Bin * as.factor(Bpath),
  family = negative.binomial(10),
  offset = offset_values
)
summary(mod_w_offset)

# Summary with and without
summary(mod_w_offset)
summary(mod)

# AIC, difference is more than 2 units
AIC(mod_w_offset)
AIC(mod)

# Compute Intercept model with offset
mod0_w_offset <- glm(Bgen9 ~ 1,
  family = negative.binomial(10),
  offset = offset_values
)

# Is Model significant than intercept model?
anova(mod0_w_offset, mod_w_offset, test = "Chisq")$Pr[2] < 0.05

# Compute R2
R2(mod)
R2(mod_w_offset)

# Now why do we see the R2 as 0.897?

# Let's first plot the model using scMaSigPro function
plotTrend(scmp.ob, "Gene9", logs = TRUE, logType = "log")

# Now lets see what is modeled
showPoly(scmp.ob)
# "beta0 + beta1*Path2vsPath1 + beta2*scmp_binned_pseudotime + beta3*scmp_binned_pseudotimexPath2"

# The data is stored in the predictor_matrix slot
predictor_matrix <- scmp.ob@Design@predictor_matrix %>% as.matrix()
print(predictor_matrix)

# Now let's check which of the terms are significant for gene-9
gene_9_sol <- showSol(scmp.ob) %>%
  rownames_to_column("gene_id") %>%
  rename(p_value = "p-value", r_squared = "R-squared") %>%
  filter(gene_id == "Gene9")
print(gene_9_sol)

# This means that the p.valor_scmp_binned_pseudotime and p.valor_scmp_binned_pseudotimexPath2 are significant
# Let's subset the predictor matrix to only include these two columns
predictor_matrix_subset <- predictor_matrix[, c("scmp_binned_pseudotime", "scmp_binned_pseudotimexPath2")]

# Since we are using the offset, we will also add the offset here
predictor_matrix_subset <- cbind(predictor_matrix_subset, offset_values) %>% as.data.frame()
print(predictor_matrix_subset)

# Now let's add the response
predictor_matrix_subset <- cbind(predictor_matrix_subset, Bgen9) %>% as.data.frame()

# Model with offset
scmp_mod_g9 <- glm(
  formula = Bgen9 ~ scmp_binned_pseudotime + scmp_binned_pseudotimexPath2,
  family = negative.binomial(10),
  data = predictor_matrix_subset, offset = predictor_matrix_subset$offset_values,
  control = list(maxit = 100), epsilon = 1e-08
)

R2(scmp_mod_g9)
## And we get the value as 0.897 which is also shown in the plot.

## Now Let's once again see the summary of "mod_w_offset"
summary(mod_w_offset)

## The "as.factor(Bpath)2" is not significant, which means we can just remove it from the model
mod_w_offset_optimized <- glm(Bgen9 ~ Bin + Bin:as.factor(Bpath),
  family = negative.binomial(10),
  offset = offset_values
)

# Plot Summary
summary(mod_w_offset_optimized)

R2(mod_w_offset_optimized)

# We have the same R2: 0.896407

## Another thing is the plot itself shows too much variance for Pseudotime 3 and 4 for path1, cosidering we have such an high r2 of 0.897
## The reason for this is the limits of the y-axis, the plot is zoomed in.
## A simple approach is to just zoom out the y-axis.

# Let's plot the data again
plotTrend(scmp.ob, "Gene9", logs = TRUE, logType = "log")

# Since this a ggplot2 object we can scale/expand the limits easily,
# let's rescale the limits of the y-axis
plotTrend(scmp.ob, "Gene9", logs = TRUE, logType = "log") +
  scale_y_continuous(
    limits = c(5, 8),
    breaks = seq(5, 8, 0.1)
  )

## Now we see the difference at the y-axis is not that much, and the plot is more informative
## Now in this perspective, the R2 of 0.897 makes more sense
## The variability itself is not more than 1 unit. If we look at Path1 at Psesudotime 3 (y = 6.1) and 4 (7.3), the model line is at 6.6-6.7.
