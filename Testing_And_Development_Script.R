# This document tempreory and is inteneded for testing and debugging
# This will be removed from the final package

# Use while debugging
#load_package()

# Set seed
set.seed(123)

save(sim.sce, file = "data/Sim2Path.RData", compress = "xz")


# Load ScMaSigpro
library(scMaSigPro)

# Step-1: Load a dataset for testing
data("Sim2Path", package = "scMaSigPro")

# Step-2: Convert to ScMaSigpro Object
scmp <- as_scmp(object = sim.sce,
                from = "sce", 
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

# Step-9: Plot the significant Genes
PlotGroups(data = scmp@compress.sce@assays@data@listData$bulk.counts[1, ],
           edesign = scmp@edesign@edesign,
           show.lines = T,
           show.fit = T,
           dis = scmp@edesign@dis,
           groups.vector = scmp@scPVector@groups.vector,
           summary.mode = "median"
           )

sc.PlotGroups(scmpObj = scmp,
              feature_id = "Gene10")
