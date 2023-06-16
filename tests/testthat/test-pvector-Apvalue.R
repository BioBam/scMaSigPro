library(testthat)
library(scMaSigPro)

test_that("Check the Adjusted p-value vector after the P-Vector Tests; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "pvector_maSigPro.RDS", package = "scMaSigPro")
  p.ob <- readRDS(data_file)

  # Step-2: Convert to a named vector
  aPvalue_vector_maSigPro <- p.ob$p.adjusted
  names(aPvalue_vector_maSigPro) <- rownames(p.ob$p.vector)

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_pVector().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract the Adjusted P-Value vectors
  aPvalue_vector_scMaSigPro <- scmp.obj@pVector@adj.p.value

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(length(aPvalue_vector_scMaSigPro), length(aPvalue_vector_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(aPvalue_vector_scMaSigPro, aPvalue_vector_maSigPro)
})
