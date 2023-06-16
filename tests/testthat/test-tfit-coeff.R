library(testthat)
library(scMaSigPro)

test_that("Check the Values for model pvalues; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "tfit_maSigPro.RDS", package = "scMaSigPro")
  tfit <- readRDS(data_file)

  # Step-2: Extract the Covariate Frame
  sol_frame_maSigPro <- tfit$sol

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_tFit().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract the path vectors
  sol_frame_scMaSigPro <-  as.data.frame(showSol(scmp.obj, view = F))

  ######################################################
  # Transfer column and row names for testing
  colnames(sol_frame_scMaSigPro) <- colnames(sol_frame_maSigPro)
  ######################################################

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(dim(sol_frame_scMaSigPro), dim(sol_frame_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(sol_frame_scMaSigPro, sol_frame_maSigPro)
})
