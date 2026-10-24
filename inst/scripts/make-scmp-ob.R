# make-scmp-ob.R
#
# Provenance script for the bundled dataset `scmp.ob` (see ?scmp.ob).
# This script is NOT run during package build or check; it documents how the
# object in data/scmp.ob.RData was produced so it can be regenerated.
# Requires only scMaSigPro. Reproduces the Quick Start vignette workflow.
# Run from the package root: Rscript inst/scripts/make-scmp-ob.R
# library() loads the installed scMaSigPro, which must be the version being
# released. To use the source tree instead, run
#   Rscript -e 'pkgload::load_all(); source("inst/scripts/make-scmp-ob.R")'
# (library() is then a no-op). The object records that version in its
# `version` slot (printed at the end); regenerate it whenever the class
# changes, so that updateObject(scmp.ob) has nothing to do.

library(scMaSigPro)

set.seed(123)
data("splat.sim", package = "scMaSigPro")

# Convert the annotated SingleCellExperiment to a ScMaSigPro object
scmp.ob <- as_scmp(
  object = splat.sim, from = "sce",
  verbose = FALSE,
  additional_params = list(
    labels_exist = TRUE,
    exist_ptime_col = "Step",
    exist_path_col = "Group"
  )
)

# Pseudo-bulk along pseudotime, polynomial degree 1, fit, select, cluster
scmp.ob <- sc.squeeze(scmp.ob)
scmp.ob <- sc.set.poly(scmp.ob, poly_degree = 1)
scmp.ob <- sc.p.vector(scmp.ob,
  offset = TRUE, p_value = 0.05, verbose = FALSE
)
scmp.ob <- sc.t.fit(scmp.ob, verbose = FALSE)
scmp.ob <- sc.filter(
  scmpObj = scmp.ob,
  rsq = 0.7,
  vars = "groups",
  intercept = "dummy",
  includeInflu = TRUE
)
scmp.ob <- sc.cluster.trend(scmp.ob, geneSet = "union", k = 4)

# Save and compress
save(scmp.ob, file = "data/scmp.ob.RData")
tools::resaveRdaFiles("data/scmp.ob.RData")
cat("scmp.ob saved by scMaSigPro", scmp.ob@version, "\n")
