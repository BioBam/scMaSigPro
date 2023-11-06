# This document is temporary and intended for testing and debugging
# This will be removed from the final package

# Set seed
set.seed(123)

# Load ScMaSigpro
library(scMaSigPro)

# Step-1: Load a dataset for testing
data("Sim2Path", package = "scMaSigPro")

# Step-2: Convert to ScMaSigpro Object
scmp.sce <- as_scmp(object = sim.sce, from = "sce",
                align_pseudotime = T,
                verbose = F,
                additional_params = list(labels_exist = TRUE,
                                         existing_pseudotime_colname = "Step",
                                         existing_path_colname = "Group")
                )

# Step-3: Pseudo-Bulk
scmp.sce <- squeeze(scmp.sce,
                         split_bins = T,
                         prune_bins = T,
                    drop_trails = T,
                         additional_params = list(use_unique_time_points = T),
                         verbose = F,
                    fill_gaps = T,
                         drop.fac = 1)

showParams(scmp.sce, view = F, return = T)

# Validation Plots
sc.plot.bins.tile(scmp.sce)
sc.plot.bins.bar(scmp.sce)

# Step-4: Make Design-Matrix
scmp <- sc.make.design.matrix(scmp.sce, poly_degree = 2)

# Step-5: Run P-vector
scmp <- sc.p.vector(scmp, parallel = T, family = gaussian())

# Step-6: Run T.fit
scmp <- sc.T.fit(scmp, parallel = T, verbose = T)

# Step-7: Select with R2 
scmp <- sc.get.siggenes(scmpObj = scmp,
                        vars = "all",
                        significant.intercept = "dummy")

nrow(showSol(scmp, view = F, return = T))
# Step-8: Plot Gene Trends
sc.PlotGroups(scmpObj = scmp,
              feature_id = "Gene435", smoothness = 0.1,
              logs = F,
              logType = "log10")


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
load("../scMaSigPro_Supp/Analysis_Public_Data/data/rep2/rep2_processed.RData")

# Convert to scmp object
scmp.cds.test.again <- as_scmp(cds,
                    "cds",
                    interactive = T,
                    verbose = F,
                    annotation_colname = "predicted.celltype.l2",
                    align_pseudotime = T)

# Bin
scmp.cds.test <- squeeze(scmp.cds.test.again,
                         split_bins = F,
                         prune_bins = F,
                         drop_trails = F,
                        additional_params = list(use_unique_time_points = T),
                        verbose = F,
                        drop.fac = 0.7)

sc.plot.bins.tile(scmpObj = scmp.cds.test)

# Validation Plots
sc.plot.bins.bar(scmpObj = scmp.cds.test)

# Step-4: Make Design-Matrix
scmp.cds.test <- sc.make.design.matrix(scmp.cds.test, poly_degree = 2)

# Step-5: Run P-vector
scmp.cds.test <- sc.p.vector(scmp.cds.test, parallel = T,
                             min.obs = 0, offset = T)

# Step-6: RunT-step
scmp.cds.test <- sc.T.fit(scmpObj = scmp.cds.test,
                     parallel = T,offset = T)

# Step-7: Get significant genes
scmp.cds.test <- sc.get.siggenes(scmpObj = scmp.cds.test, vars = "groups")


# View Sig genes
View(showSol(scmp.cds.test, return = T, view = F))

# Plots
sc.PlotGroups(scmpObj = scmp.cds.test,
              feature_id = "HUWE1",
              smoothness = 2,
              logs = T,
              logType = "log")

sc.PlotGroups(scmpObj = scmp.cds.test,
              feature_id = "CDK14",
              smoothness = 0.01,
              logs = T,
              logType = "log")
stop()
