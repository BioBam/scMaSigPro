suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Check-'fit$p.adjusted' Reproducibility: Match, Dimension, Name and Value", {
  # Step-1: Load Data
  data("data.abiotic")
  data("edesign.abiotic")

  # Create MaSigPro Design
  design_2 <- make.design.matrix(edesign = edesign.abiotic, degree = 2)
  design_3 <- make.design.matrix(edesign = edesign.abiotic, degree = 3)
  design_4 <- make.design.matrix(edesign = edesign.abiotic, degree = 4)

  # Run-p-vector
  gc <- capture.output(
    fit2 <- p.vector(data.abiotic, design_2, Q = 0.05, MT.adjust = "BH", min.obs = 20)
  )
  gc <- capture.output(
    fit3 <- p.vector(data.abiotic, design_3, Q = 0.05, MT.adjust = "BH", min.obs = 20)
  )
  gc <- capture.output(
    fit4 <- p.vector(data.abiotic, design_4, Q = 0.05, MT.adjust = "BH", min.obs = 20)
  )

  # Convert the Design for scMaSigPro
  add_col <- function(x) {
    return(names(x[x == 1]))
  }
  edesign.abiotic.in <- as.data.frame(edesign.abiotic[, 1, drop = F])
  Group <- c(unlist(apply(edesign.abiotic[, c(3:6)], 1, add_col, simplify = F), use.names = F))
  edesign.abiotic.in$Group <- Group

  # Create Design with scMaSigPro
  test.scmp <- create_scmpObj(
    counts = data.abiotic,
    cell_data = edesign.abiotic.in,
    pseudotime_colname = "Time",
    path_colname = "Group",
    use_as_bin = T
  )

  # Run Make Design
  test.scmp.2 <- sc.make.design.matrix(test.scmp, poly_degree = 2)
  test.scmp.3 <- sc.make.design.matrix(test.scmp, poly_degree = 3)
  test.scmp.4 <- sc.make.design.matrix(test.scmp, poly_degree = 4)

  # Run sc pvectory
  test.scmp.2 <- sc.p.vector(test.scmp.2,
    min.obs = 20,
    offset = F, parallel = T,
    family = gaussian(), verbose = F
  )

  test.scmp.3 <- sc.p.vector(test.scmp.3,
    min.obs = 20,
    offset = F, parallel = T,
    family = gaussian(), verbose = F
  )

  test.scmp.4 <- sc.p.vector(test.scmp.4,
    min.obs = 20,
    offset = F, parallel = T,
    family = gaussian(), verbose = F
  )

  # Matric Convert
  names(test.scmp.2@scPVector@FDR) <- NULL
  names(test.scmp.3@scPVector@FDR) <- NULL
  names(test.scmp.4@scPVector@FDR) <- NULL

  # Check
  # Poly-order-2
  expect_equal(
    expected = fit2$FDR,
    object = test.scmp.2@scPVector@FDR
  )
  # Poly-order-3
  expect_equal(
    expected = fit3$FDR,
    object = test.scmp.3@scPVector@FDR
  )
  # Poly-order-4
  expect_equal(
    expected = fit4$FDR,
    object = test.scmp.4@scPVector@FDR
  )
})
