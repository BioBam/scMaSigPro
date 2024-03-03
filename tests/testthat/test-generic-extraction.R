suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))

test_that("Check Extraction Generics", {
  # Step-1: Load Data
  data("splat.sim", package = "scMaSigPro")

  # Step-2: Create scMaSigPro object
  scmp_ob <- as_scmp(
    object = splat.sim, from = "sce",
    align_pseudotime = FALSE,
    verbose = FALSE,
    additional_params = list(
      labels_exist = TRUE,
      exist_ptime_col = "Step",
      exist_path_col = "Group"
    )
  )

  # Step-3: Perform Pseudobulking
  scmp_ob <- sc.squeeze(scmp_ob)

  # Check
  expect_equal(
    expected = as.matrix(scmp_ob@Sparse@assays@data@listData$counts),
    object = eSparse(scmp_ob)
  )

  expect_equal(
    expected = as.matrix(scmp_ob@Dense@assays@data@listData$bulk.counts),
    object = eDense(scmp_ob)
  )
})
