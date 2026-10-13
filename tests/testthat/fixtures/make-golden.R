# Regenerates tests/testthat/fixtures/golden-v0.0.4.rds, the frozen outputs
# that test-golden.R compares against.
#
# Run from the package root:
#     Rscript tests/testthat/fixtures/make-golden.R
#
# Regenerate only when a change is MEANT to alter results. Commit the new file
# in the same commit as the code change, summarise what moved in the commit
# message, and add a NEWS entry under "Changes affecting results".
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-fixtures.R")

set.seed(2024)
golden <- build_golden()

out <- "tests/testthat/fixtures/golden-v0.0.4.rds"
saveRDS(golden, out, compress = "xz")
cat("wrote", out, "(", round(file.size(out) / 1024), "KB )\n")
