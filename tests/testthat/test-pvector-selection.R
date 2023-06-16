library(testthat)
library(scMaSigPro)

test_that("Check the result from the P-vector selection at the same a-p-value; match dimension and value", {
  # Step-1: Load object from /extdata
  data_file <- system.file("extdata", "pvector_maSigPro.RDS", package = "scMaSigPro")
  p.ob <- readRDS(data_file)

  # Step-2: Extract the Selected Genes
  sel_pvector_maSigPro <- as.matrix(p.ob$SELEC)

  # Step-4: Load the object of scMaSigPro object after create_scMaSigPro_obj()
  data_file <- system.file("extdata", "sc_pVector().RDS", package = "scMaSigPro")
  scmp.obj <- readRDS(data_file)

  # Step-5: Extract significance threshold
  sig.thresh <- p.ob$Q

  # Extract the names of the genes
  sig.gene.names <- names(scmp.obj@pVector@adj.p.value[scmp.obj@pVector@adj.p.value <= scmp.obj@parameters@p.vector.sig])

  # Extract the frame
  sel_pvector_scMaSigPro <- as.matrix(scmp.obj@assays@data@listData$counts[rownames(scmp.obj@assays@data@listData$counts) %in% sig.gene.names, , drop = F])

  # Step-6: First, check that the matrices have the same dimensions
  expect_equal(dim(sel_pvector_scMaSigPro), dim(sel_pvector_maSigPro))

  # Step-7: Next, check that all elements in the matrices are the same
  expect_equal(sel_pvector_scMaSigPro, sel_pvector_maSigPro)
})
