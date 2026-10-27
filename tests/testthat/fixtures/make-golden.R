# Regenerates tests/testthat/fixtures/golden-v0.0.4.rds, the frozen outputs
# that test-golden.R compares against.
#
# Run from the package root:
#     Rscript tests/testthat/fixtures/make-golden.R
#
# Regenerate only when a change is MEANT to alter results. Commit the new file
# in the same commit as the code change, summarise what moved in the commit
# message, and add a NEWS entry under "Changes affecting results".
# testthat runs every test file with LC_COLLATE=C; generate under the same
# collation so sort()/order() on labels cannot differ.
invisible(Sys.setlocale("LC_COLLATE", "C"))
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-fixtures.R")

# build_golden() must be deterministic: test-golden.R sets no seed.
golden <- build_golden()

out <- "tests/testthat/fixtures/golden-v0.0.4.rds"
saveRDS(golden, out, compress = "xz")
cat("wrote", out, "(", round(file.size(out) / 1024), "KB )\n")
