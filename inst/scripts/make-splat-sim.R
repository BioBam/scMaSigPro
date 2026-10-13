# make-splat-sim.R
#
# Provenance script for the bundled dataset `splat.sim` (see ?splat.sim).
# This script is NOT run during package build or check; it documents how the
# object in data/splat.sim.RData was produced so it can be regenerated.
# Requires: splatter, scran, scuttle, scater, SingleCellExperiment.
# Run from the package root: Rscript inst/scripts/make-splat-sim.R

# Load Required Packages
suppressPackageStartupMessages(library(splatter))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(SingleCellExperiment))

set.seed(123)

# Simulate
splat.sim <- splatSimulatePaths(
   params = newSplatParams(
       batchCells = 200, nGenes = 100),
   group.prob = c(0.5, 0.5),
   path.nSteps = c(100, 100),
   de.prob = 0.3, de.facLoc = 0.2,
   path.from = c(0, 0), # Bifurcation
   verbose = FALSE)

# Normalize
 splat.sim <- logNormCounts(splat.sim, assay.type = "counts")

# Reduce Dimensions
 splat.sim <- runPCA(splat.sim, exprs_values = "logcounts", ncomponents = 2)

# Visulize Steps and Groups
 plotPCA(splat.sim, colour_by = "Step")
 plotPCA(splat.sim, colour_by = "Group")

# Create sce and transfer data
sce <- SingleCellExperiment(list(counts = splat.sim@assays@data@listData$counts))
sce@colData <- splat.sim@colData
rowData(sce) <- rowData(splat.sim)
reducedDims(sce) <- reducedDims(splat.sim)
splat.sim <- sce

# Save
 save(splat.sim, file = "data/splat.sim.RData")

# Compress
tools::resaveRdaFiles("data/splat.sim.RData")
