# This document tempreory and is inteneded for testing and debugging
# This will be removed from the final package

# Set seed
set.seed(123)

# Load ScMaSigpro
library(scMaSigPro)

# Step-1: Load a dataset for testing
data("Sim2Path", package = "scMaSigPro")

# Step-2: Convert to ScMaSigpro Object
scmp <- as_scmp(object = sim.sce,
                from = "sce",
                additional_params = list(
                    overwrite_labels = TRUE,
                    existing_pseudotime_colname = "Step",
                    existing_path_colname = "Group"
                ))

# Step-3-A: Perform Binning
scmp <- entropy_discretize(scmp)

# Step-3-B: Create Pseudo-Bulk Cell Metadata
scmp <- make.pseudobulk.design(scmp)

# Step-3-C: Pseudo-bulk the counts
scmp <- make.pseudobulk.counts(scmp)

# Step-4: Run Pvector

