suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))

test_that("Check Pseudo-bulking with manual bulking", {
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

  # Extract data for a random gene
  random_gene <- sample(nrow(eSparse(scmp_ob)), 1)
  slot_gen9 <- eSparse(scmp_ob)[random_gene, ]
  slot_Bgen9 <- eDense(scmp_ob)[random_gene, ]

  # Extract Pseudotime and the binned Pseudotime
  Pseudotime <- cSparse(scmp_ob)$Pseudotime
  names(Pseudotime) <- rownames(cSparse(scmp_ob))
  BinnedPseudotime <- cSparse(scmp_ob)$scmp_binned_pseudotime
  names(BinnedPseudotime) <- rownames(cSparse(scmp_ob))

  # Extract bin information
  Bin <- cDense(scmp_ob)[[scmp_ob@Parameters@bin_ptime_col]]
  names(Bin) <- rownames(cDense(scmp_ob))
  Bpath <- rep(c(1:2), each = length(Bin) / 2)

  # Check total sum
  expect_identical(
    expected = as.numeric(sum(slot_gen9)),
    object = as.numeric(sum(slot_Bgen9))
  )

  # Create vectors
  sum1_new <- sum2_new <- NULL

  # Compute pseudo-bulk
  for (i in 1:length(Bin) / 2) {
    cell_vector <- names(BinnedPseudotime[BinnedPseudotime == i])
    cell_count_vector <- slot_gen9[names(slot_gen9) %in% cell_vector]
    sum1_new <- c(sum1_new, sum(cell_count_vector))

    bin_vector <- names(Bin[Bin == i])
    bin_count_vector <- slot_Bgen9[names(Bin) %in% bin_vector]
    sum2_new <- c(sum2_new, sum(bin_count_vector))
  }

  # Check
  expect_identical(
    object = as.numeric(sum1_new),
    expected = as.numeric(sum2_new)
  )
})
