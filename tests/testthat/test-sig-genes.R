library(testthat)
library(scMaSigPro)

test_that("Check Selections of significant genes; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "sigs_maSigPro.RDS", package = "scMaSigPro")
  sigs <- readRDS(data_file)

  # Step-2: Extract the Covariate Frame
  summary_frame_maSigPro <- sigs$summary

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "scmp_sigGenes_sc.sel.features().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract the path vectors
  summary_frame_scMaSigPro <- as.data.frame(scmp.obj@sigGenes@sel.genes.frame)

  ######################################################
  # Transfer column shift columns
  summary_frame_scMaSigPro <- summary_frame_scMaSigPro[, c(2,1),drop = F]
  colnames(summary_frame_scMaSigPro) <- colnames(summary_frame_maSigPro)
  ######################################################

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(dim(summary_frame_scMaSigPro), dim(summary_frame_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(summary_frame_scMaSigPro, summary_frame_maSigPro)
})
