suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Check-'design$dis' Reproducibility", {
  # Step-1: Load Data
  data("data.abiotic")
  data("edesign.abiotic")

  # Step-2: Set-up data for scMaSigPro
  count <- as.matrix(data.abiotic)
  cell_metadata <- as.data.frame(edesign.abiotic)

  # Step-2.1: Add group column
  cell_metadata$Group <- apply(cell_metadata[, c(3:6)], 1, FUN = function(x) {
    return(names(x[x == 1]))
  })

  # Step-2.2: Remove Binary Columns
  cell_metadata <- cell_metadata[, !(colnames(cell_metadata) %in%
    c("Control", "Cold", "Heat", "Salt")),
  drop = FALSE
  ]

  # Step-3: Create scmp Object
  test.scmp <- create.scmp(
    counts = count,
    cell_data = cell_metadata,
    pseudotime_colname = "Time",
    path_colname = "Group",
    use_as_bin = T
  )

  # Step-4: Create MaSigPro Design
  design_2 <- make.design.matrix(edesign = edesign.abiotic, degree = 2)
  design_3 <- make.design.matrix(edesign = edesign.abiotic, degree = 3)
  design_4 <- make.design.matrix(edesign = edesign.abiotic, degree = 4)

  # Step-5: Set polynomial and make design
  test.scmp.2 <- sc.set.poly(test.scmp, poly_degree = 2)
  test.scmp.3 <- sc.set.poly(test.scmp, poly_degree = 3)
  test.scmp.4 <- sc.set.poly(test.scmp, poly_degree = 4)

  # Check-dis
  # Poly-order-2
  expect_identical(
    expected = design_2$dis,
    object = test.scmp.2@design@predictor_matrix
  )
  # Poly-order-3
  expect_identical(
    expected = design_3$dis,
    object = test.scmp.3@design@predictor_matrix
  )
  # Poly-order-4
  expect_identical(
    expected = design_4$dis,
    object = test.scmp.4@design@predictor_matrix
  )
})
