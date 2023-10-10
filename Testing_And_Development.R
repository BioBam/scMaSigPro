##########################################################
## Author: Priyansh Srivastava ###########################
## Email: spriyansh29@gmail.com ##########################
## Script: scMaSigPro Application ########################
##########################################################

set.seed(007)

# Load scMaSigPro
library(scMaSigPro)
library(ggplot2)

# Prefix
prefixIn <- "../scMaSigPro_Supp/benchmarks/01_Sparsity/data/simulated/sce/"

# Load CDS object
load(paste0(prefixIn, "sparsity_30.RData"))

# Monocl3 3 object
sim.sce

# Convert single-cell experiment object
scmp.obj <- as_scmp(sim.sce, from = "sce")

# All binning in one go
scmp.obj <- squeeze(
    scmp.ob = scmp.obj,
    time.col = "Step",
    path.col = "Group",
    method = "Sturges",
    drop.fac = 1,
    verbose = T,
    cluster.count.by = "mean"
)

# Make Design
scmp.obj <- sc.make.design.matrix(scmp.obj,
                                  degree = 2,
                                  time.col = "binnedTime",
                                  path.col = "path"
)

# Run p-vector
scmp.obj <- sc.p.vector(
        scmpObj = scmp.obj, verbose = T, min.obs = 6,
        counts = T, theta = 1,
        offset = T, epsilon = 0.00001
        )

# Run-Step-2
scmp.obj <- sc.T.fit(
        data = scmp.obj, verbose = T,
        step.method = "backward",
        family = scmp.obj@scPVector@family,
        offset = T
)

# Extract the genes
scmp.obj <- sc.get.siggenes(scmpObj = scmp.obj, rsq = 0.7, vars = "groups")

# Plot the data
sc.PlotGroups(scmpObj = scmp.obj,feature_id = "Gene45", dis = scmp.obj@scTFit@dis,
                       edesign =  scmp.obj@scTFit@edesign,
                       groups.vector = scmp.obj@scTFit@groups.vector)

sc.plot.bins(scmpObj = scmp.obj)


sc.path.intersection(scmpObj = scmp.obj)



saveRDS(scmp.obj, "../scMaSigPro_Supp/Testing_And_Development/scmp.obj.latest.RDS")
scmp.obj <- readRDS("../scMaSigPro_Supp/Testing_And_Development/scmp.obj.latest.RDS")

