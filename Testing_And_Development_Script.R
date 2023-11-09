# This document is temporary and intended for testing and debugging
# This will be removed from the final package

# Set seed
set.seed(123)

# Load ScMaSigpro
library(scMaSigPro)

# Step-1: Load a dataset for testing
## This dataset will be avaible as the part of the package
data("Sim2Path", package = "scMaSigPro")

# Step-2: Convert to ScMaSigpro Object
## Here we will convert the SCE object to scMaSigPro object
scmp.sce <- as_scmp(
  object = sim.sce, from = "sce",
  align_pseudotime = T,
  verbose = F,
  additional_params = list(
    labels_exist = TRUE,
    existing_pseudotime_colname = "Step",
    existing_path_colname = "Group"
  )
)

# Step-3: Pseudo-Bulk
# This is the main stp
scmp.sce <- squeeze(
  scmpObject = scmp.sce,
  bin_method = "Sturges",
  drop.fac = 1,
  verbose = F,
  cluster_count_by = "sum",
  split_bins = FALSE,
  prune_bins = F,
  drop_trails = F,
  fill_gaps = F
)

# Validation Plots
# sc.plot.bins.tile(scmp.sce)
# sc.plot.bins.bar(scmp.sce)

# Step-4: Make Design-Matrix
scmp.sce <- sc.make.design.matrix(scmp.sce, poly_degree = 2)

# Step-5: Run P-vector
# offset_F_UseWeights_F_UseInverseWeights_F_UseBinWeightAsOffset_T
scmp.sce <- sc.p.vector(scmp.sce,
  parallel = F, useWeights = F,
  offset = F, useInverseWeights = F, min.obs = 1,
  logOffset = F, globalTheta = F
)

# Step-6: Run T.fit
scmp.sce <- sc.T.fit(scmp.sce, parallel = T, verbose = T)


showParams(scmp.sce, return = T, view = F)

# Step-7: Select with R2
scmp.sce <- sc.get.siggenes(
  scmpObj = scmp.sce,
  vars = "all",
  significant.intercept = "dummy"
)
showParams(scmp.sce, return = T, view = F)

nrow(showSol(scmp.sce, view = F, return = T, influ = F))
# Step-8: Plot Gene Trends
sc.PlotGroups(
  scmpObj = scmp.sce,
  feature_id = "Gene121", smoothness = 0.1,
  logs = F,
  logType = "log2"
)


# Developing Methods for monocle3
# test real Data
library(shiny)
library(plotly)
library(assertthat)
library(tidyverse)
library(SingleCellExperiment)
library(entropy)
library(igraph)

# Step-1: Load data and Monocle3 like function
load("../scMaSigPro_Supp/Analysis_Public_Data/data/rep3/rep3_processed.RData")

# Convert to scmp object
scmp.cds.test <- as_scmp(cds,
  "cds",
  interactive = T,
  verbose = F,
  annotation_colname = "predicted.celltype.l2",
  align_pseudotime = F
)

# Bin
scmp.cds.test <- squeeze(scmp.cds.test,
  split_bins = T,
  prune_bins = T,
  drop_trails = T,
  additional_params = list(use_unique_time_points = T),
  verbose = F,
  drop.fac = 1
)

sc.plot.bins.tile(scmpObj = scmp.cds.test)

# Validation Plots
sc.plot.bins.bar(scmpObj = scmp.cds.test)

# Step-4: Make Design-Matrix
scmp.cds.test <- sc.make.design.matrix(scmp.cds.test, poly_degree = 2)

# Step-5: Run P-vector
scmp.cds.test <- sc.p.vector(scmp.cds.test,
  parallel = T,
  min.obs = 0, offset = T
)

# Step-6: RunT-step
scmp.cds.test <- sc.T.fit(
  scmpObj = scmp.cds.test,
  parallel = T, offset = T
)

# Step-7: Get significant genes
scmp.cds.test <- sc.get.siggenes(scmpObj = scmp.cds.test, vars = "groups")


# View Sig genes
View(showSol(scmp.cds.test, return = T, view = F, influ = F))

# Plots
sc.PlotGroups(
  scmpObj = scmp.cds.test,
  feature_id = "SORBS3",
  smoothness = 2,
  logs = T,
  logType = "log"
)

sc.PlotGroups(
  scmpObj = scmp.cds.test,
  feature_id = "CDK14",
  smoothness = 0.01,
  logs = T,
  logType = "log"
)
stop()



### Testing for Parameterization of Negative Binomial
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(glmnet))


# Get gene metadata
gene.metadata <- SingleCellExperiment::rowData(scmp.sce@sce) %>% as.data.frame()

# Get Matrix
data <- scmp.sce@scPVector@dis

# Get response
y_matrix <- as.matrix(scmp.sce@compress.sce@assays@data@listData$bulk.counts)

# Select certain kinds of genes
y_sel_df <- data.frame(
  geneName = c("Gene476", "Gene1", "Gene391", "Gene261", "Gene85", "Gene138", "Gene99"),
  outlier = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE),
  mean = c("High", "High", "High", "Low", "Low", "Low", "Low"),
  pathChange = c(
    "DownPath2", "NoChange", "UpPath2",
    "UpPath1", "NoChange", "DownPath2",
    "NoChange"
  )
)

# Extract genes
y_sel_matrix <- y_matrix[y_sel_df$geneName, ]
offsets_data <- estimateSizeFactorsForMatrix(y_sel_matrix)

library(MASS)

# Function to fit a model and return the model object
fit_nb_model <- function(gene_expression, covars_df, offsets) {
  # The gene_expression vector contains the counts for a single gene across different bins
  covars_df[["y"]] <- gene_expression

  # Fit the negative binomial model
  model <- glm.nb(y ~ ., data = covars_df)

  return(model)
}

# Initialize a list to store the models for each gene
model.list <- list()

# Loop over the rows of y_sel to fit a model for each gene
for (gene_i in rownames(y_sel_matrix)) {
  # Extract the gene expression data for the current gene
  gene_expression <- as.numeric(y_sel_matrix[gene_i, ])

  # Fit the model and store it in the list
  model.list[[gene_i]] <- fit_nb_model(gene_expression, data, offsets = log(offsets_data))
}

# Now you have a list 'model.list' that contains the fitted models for each gene
# You can extract the estimated theta for each gene using sapply:
theta.list <- sapply(model.list, function(model) model$theta)

# The 'theta.list' object now contains the theta values for each gene
