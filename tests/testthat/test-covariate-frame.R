library(testthat)
library(scMaSigPro)

test_that("Check the covariate dataframe; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "tfit_maSigPro.RDS", package = "scMaSigPro")
  tfit <- readRDS(data_file)

  # Step-2: Extract the Covariate Frame
  covariate_frame_maSigPro <- tfit$dis

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_makeDesign().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract the path vectors
  covariate_frame_scMaSigPro <- as.data.frame(scmp.obj@covariate@covariate)

  ######################################################
  # Transfer column and row names for testing
  colnames(covariate_frame_scMaSigPro) <- colnames(covariate_frame_maSigPro)
  rownames(covariate_frame_scMaSigPro) <- rownames(covariate_frame_maSigPro)
  ######################################################

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(dim(covariate_frame_scMaSigPro), dim(covariate_frame_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(covariate_frame_scMaSigPro, covariate_frame_maSigPro)
})
