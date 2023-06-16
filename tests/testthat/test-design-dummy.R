library(testthat)
library(scMaSigPro)

test_that("Dummy variables match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "tfit_maSigPro.RDS", package = "scMaSigPro")
  tfit <- readRDS(data_file)

  # Step-2: Extract the Design file
  design.matrix <- tfit$edesign

  # Step-3: Get the vectors for the path
  path1_vec_maSigPro <- design.matrix$Path1
  path2_vec_maSigPro <- design.matrix$Path2

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_makeDesign().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-6: Extract the path vectors
  path1_vec_scMaSigPro <- scmp.obj@covariate@factor.design[, 1]
  path2_vec_scMaSigPro <- scmp.obj@covariate@factor.design[, 2]

  # Step-7: First, check that the matrices have the same dimensions
  expect_equal(length(path1_vec_maSigPro), length(path1_vec_scMaSigPro))
  expect_equal(length(path2_vec_maSigPro), length(path2_vec_scMaSigPro))

  # Step-8: Next, check that all elements in the matrices are the same
  expect_equal(path1_vec_maSigPro, path1_vec_scMaSigPro)
  expect_equal(path2_vec_maSigPro, path2_vec_scMaSigPro)
})
