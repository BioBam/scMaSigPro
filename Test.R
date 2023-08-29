##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: scMaSigPro Application ########################
##########################################################

set.seed(007)

# Call the required libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(scMaSigPro))


# Prefix
prefixIn <- "../scMaSigPro_Supp/benchmarks/11_RealDataSmall/data/results/"
prefixOut <- "../scMaSigPro_Supp/benchmarks/11_RealDataSmall/data/output/"

# Load CDS object
load(paste0(prefixIn, "monocle3_inferred_pseudotime.RData"))

# Convert the ScMaSigPro Object
scmp.obj <- as_scmp(cds, from = "cell_data_set")

# Plot the Paths
plot_cells(scmp.obj@sce, color_cells_by = "Path")

# Subset
scmp.obj_test <- scMaSigPro::selectPath(obj = scmp.obj, sel.path = c("Path1", "Path2"),
                            balance_paths = T, pathCol = "Path", pTimeCol = "Pseudotime",
                            plot_paths = F, verbose = F)

# Plot the Paths
plot_cells(scmp.obj@sce, color_cells_by = "Path")

# Compress
scmp.obj <- squeeze(
    scmp.ob = scmp.obj,
    time.col = "Pseudotime",
    path.col = "Path",
    method = "Sturges",
    drop.fac = 0.7,
    verbose = T,
    cluster.count.by = "sum",
)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  degree = 2,
                                  time.col = "binnedTime",
                                  path.col = "path"
)

# Run p-vector
scmp.obj <- sc.p.vector(
    scmpObj = scmp.obj, verbose = T, min.obs = 10,
    counts = T, theta = 10,
    offset = T
)

# Run-Step-2
scmp.obj <- sc.T.fit(
    data = scmp.obj, verbose = T,
    step.method = "backward",
    family = scmp.obj@scPVector@family,
    offset = T
)

# Extract the genes
sig.gene <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")


sig.gene$summary[, 1][1]


blk.counts <- scmp.obj@compress.sce@assays@data@listData$bulk.counts

scMaSigPro::PlotGroups(data = blk.counts[sig.gene$summary[, 2][11],],
                       edesign =  scmp.obj@scPVector@edesign)