suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Check-'tstep$SOL' Reproducibility: Match, Dimension, Name and Value", {
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

  # Run-tstep
  gc <- capture.output(
    tstep2 <- T.fit(fit2, step.method = "backward", alfa = 0.05)
  )
  gc <- capture.output(
    tstep3 <- T.fit(fit3, step.method = "backward", alfa = 0.05)
  )
  gc <- capture.output(
    tstep4 <- T.fit(fit4, step.method = "backward", alfa = 0.05)
  )

  # Extract All
  sigs2 <- get.siggenes(tstep2, rsq = 0.6, vars = "all")
  sigs3 <- get.siggenes(tstep3, rsq = 0.6, vars = "all")
  sigs4 <- get.siggenes(tstep4, rsq = 0.6, vars = "all")

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

  # Run Tstep
  test.scmp.2 <- sc.T.fit(test.scmp.2,
    parallel = T, verbose = F, offset = F
  )
  test.scmp.3 <- sc.T.fit(test.scmp.3,
    parallel = T, verbose = F, offset = F
  )
  test.scmp.4 <- sc.T.fit(test.scmp.4,
    parallel = T, verbose = F, offset = F
  )

  # Extract sigs
  test.scmp.2 <- sc.get.siggenes(test.scmp.2, rsq = 0.6, vars = "all",use_Influ = T)
  test.scmp.3 <- sc.get.siggenes(test.scmp.3, rsq = 0.6, vars = "all",use_Influ = T)
  test.scmp.4 <- sc.get.siggenes(test.scmp.4, rsq = 0.6, vars = "all",use_Influ = T)

  # Check
  # Poly-order-2
  expect_equal(
    expected = sigs2$summary,
    object = test.scmp.2@sig.genes@summary
  )
  # Poly-order-3
  expect_equal(
    expected = sigs3$summary,
    object = test.scmp.3@sig.genes@summary
  )
  # Poly-order-4
  expect_equal(
    expected = sigs4$summary,
    object = test.scmp.4@sig.genes@summary
  )
})
