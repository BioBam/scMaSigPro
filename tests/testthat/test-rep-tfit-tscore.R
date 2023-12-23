suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Check-'tfit$t.score'", {
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

  # Step-6: Run P.vector
  gc <- capture_output(fit_2 <- p.vector(data.abiotic, design_2,
    Q = 0.05,
    MT.adjust = "BH", min.obs = 20
  ))
  gc <- capture_output(fit_3 <- p.vector(data.abiotic, design_3,
    Q = 0.05,
    MT.adjust = "BH", min.obs = 20
  ))
  gc <- capture_output(fit_4 <- p.vector(data.abiotic, design_4,
    Q = 0.05,
    MT.adjust = "BH", min.obs = 20
  ))

  # Step-7: Run sc.p.vector
  test.scmp.2 <- sc.p.vector(test.scmp.2,
    min.na = 20, verbose = FALSE,
    offset = FALSE, parallel = FALSE, max_it = 25,
    epsilon = 0.00001, family = gaussian()
  )
  test.scmp.3 <- sc.p.vector(test.scmp.3,
    min.na = 20, verbose = FALSE,
    offset = FALSE, parallel = FALSE, max_it = 25,
    epsilon = 0.00001, family = gaussian()
  )
  test.scmp.4 <- sc.p.vector(test.scmp.4,
    min.na = 20, verbose = FALSE,
    offset = FALSE, parallel = FALSE, max_it = 25,
    epsilon = 0.00001, family = gaussian()
  )

  # Step-8: Run t.fit
  gc <- capture_output(tstep_2 <- T.fit(fit_2, step.method = "backward", alfa = 0.05))
  gc <- capture_output(tstep_3 <- T.fit(fit_3, step.method = "backward", alfa = 0.05))
  gc <- capture_output(tstep_4 <- T.fit(fit_4, step.method = "backward", alfa = 0.05))

  # Step-9: Run sc.t.fit
  test.scmp.2 <- sc.t.fit(test.scmp.2, verbose = FALSE)
  test.scmp.3 <- sc.t.fit(test.scmp.3, verbose = FALSE)
  test.scmp.4 <- sc.t.fit(test.scmp.4, verbose = FALSE)

  # Check-sol
  # Poly-order-2
  expect_identical(as.data.frame(showTS(test.scmp.2)), expected = tstep_2$t.score)
  # Poly-order-3
  expect_identical(as.data.frame(showTS(test.scmp.3)), expected = tstep_3$t.score)
  # Poly-order-4
  expect_identical(as.data.frame(showTS(test.scmp.4)), expected = tstep_4$t.score)
})
