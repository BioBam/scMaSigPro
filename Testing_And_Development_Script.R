# This document is temporary and intended for testing and debugging
# This will be removed from the final package

# Set seed
set.seed(123)

# Load ScMaSigpro
library(scMaSigPro)
library(ggplot2)

# Step-1: Load a dataset for testing
data("Sim2Path", package = "scMaSigPro")

# Step-2: Convert to ScMaSigpro Object
scmp <- as_scmp(object = sim.sce, from = "sce",
                additional_params = list(labels_exist = TRUE,
                                         existing_pseudotime_colname = "Step",
                                         existing_path_colname = "Group")
                )

# Step-3-A: Perform Binning
scmp <- entropy_discretize(scmp,drop.fac = 1,
                           verbose = T,
                           binning = "individual",
                           additional_params = list(use_unique_time_points = TRUE))

# Step-3-B: Create Pseudo-Bulk Cell Metadata
scmp <- make.pseudobulk.design(scmp,
                               verbose= T)

# Step-3-C: Pseudo-bulk the counts
scmp <- make.pseudobulk.counts(scmp)

# Step-4: Make Design-Matrix
scmp <- sc.make.design.matrix(scmp, poly_degree = 2)

# Step-5: Run P-vector
scmp <- sc.p.vector(scmp, parallel = T, family = gaussian())

# Step-6: Run T.fit
scmp <- sc.T.fit(scmp, parallel = T, verbose = T)

# Step-7: Select with R2 
scmp <- sc.get.siggenes(scmpObj = scmp,
                        vars = "all",
                        significant.intercept = "dummy")

# Step-8: Plot Gene Trends
sc.PlotGroups(scmpObj = scmp,
              feature_id = "Gene10", smoothness = 0.1)


# Developing Methods for monocle3

# Step-1: Load data and Monocle3 like function
load("extdata/rep1_processed.RData")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(SingleCellExperiment))

# Step-2: Select the Paths
scmp.cds  <- selectPath.m3(cdsObj = cds, annotation = "predicted.celltype.l2")

# Create SCMP Object
test <-  as_scmp(cds, "cds")

scmp.cds <- entropy_discretize(scmp.cds,drop.fac = 0.4,
                           verbose = T,
                           binning = "individual",
                           additional_params = list(use_unique_time_points = TRUE))

scmp.cds <- make.pseudobulk.design(scmp.cds,
                               verbose= T)
scmp.cds <- make.pseudobulk.counts(scmp.cds)

# Step-4: Make Design-Matrix
scmp.cds <- sc.make.design.matrix(scmp.cds, poly_degree = 2)

# Step-5: Run P-vector
scmp.cds <- sc.p.vector(scmp.cds, parallel = T)

# Step-6: RunT-step
scmp.cds <- sc.T.fit(scmpObj = scmp.cds,
                     parallel = T)

# Step-7: Get significant genes
scmp.cds <- sc.get.siggenes(scmpObj = scmp.cds, vars = "groups")

# Step-8: Plot
sc.PlotGroups(scmpObj = scmp.cds,
              feature_id = "IRF8")
