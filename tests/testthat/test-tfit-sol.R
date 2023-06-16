library(testthat)
library(scMaSigPro)

test_that("Check the Values for model beta estimates; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "tfit_maSigPro.RDS", package = "scMaSigPro")
  tfit <- readRDS(data_file)

  # Step-2: Extract the Covariate Frame
  coeff_frame_maSigPro <- tfit$coefficients

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_tFit().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract the path vectors
  coeff_frame_scMaSigPro <- as.data.frame(showEstimates(scmp.obj, view = F))

  ######################################################
  # Transfer column and row names for testing
  colnames(coeff_frame_scMaSigPro) <- colnames(coeff_frame_maSigPro)
  ######################################################

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(dim(coeff_frame_scMaSigPro), dim(coeff_frame_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(coeff_frame_scMaSigPro, coeff_frame_maSigPro)
})
