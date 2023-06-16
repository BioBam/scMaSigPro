library(testthat)
library(scMaSigPro)

test_that("Check the p-value vector after the P-Vector Tests; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "pvector_maSigPro.RDS", package = "scMaSigPro")
  p.ob <- readRDS(data_file)

  # Step-2: Extract the P value
  pvalue_frame_maSigPro <- p.ob$p.vector

  # Step-3: Convert to a named vector
  pvalue_vector_maSigPro <- as.vector(pvalue_frame_maSigPro)
  names(pvalue_vector_maSigPro) <- rownames(pvalue_frame_maSigPro)

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_pVector().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract the P-Value vectors
  pvalue_vector_scMaSigPro <- scmp.obj@pVector@p.value

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(length(pvalue_vector_scMaSigPro), length(pvalue_vector_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(pvalue_vector_scMaSigPro, pvalue_vector_maSigPro)
})
