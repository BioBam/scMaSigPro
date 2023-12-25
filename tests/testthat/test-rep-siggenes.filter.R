suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Check-'sigs$clusters'", {
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
  test.scmp <- create_scmp(
    counts = count,
    cell_data = cell_metadata,
    ptime_col = "Time",
    path_col = "Group",
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
    min_na = 20, verbose = FALSE,
    offset = FALSE, parallel = FALSE, max_it = 25,
    epsilon = 0.00001, family = gaussian()
  )
  test.scmp.3 <- sc.p.vector(test.scmp.3,
    min_na = 20, verbose = FALSE,
    offset = FALSE, parallel = FALSE, max_it = 25,
    epsilon = 0.00001, family = gaussian()
  )
  test.scmp.4 <- sc.p.vector(test.scmp.4,
    min_na = 20, verbose = FALSE,
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

  # Step-10: Compute Clusters with maSigPro
  sigs.2 <- get.siggenes(tstep_2, rsq = 0.6, vars = "groups")
  sigs.3 <- get.siggenes(tstep_3, rsq = 0.6, vars = "groups")
  sigs.4 <- get.siggenes(tstep_4, rsq = 0.6, vars = "groups")

  # Step-11: Compute Clusters with scMaSigPro
  test.scmp.2 <- sc.filter(test.scmp.2, rsq = 0.6, vars = "groups")
  test.scmp.3 <- sc.filter(test.scmp.3, rsq = 0.6, vars = "groups")
  test.scmp.4 <- sc.filter(test.scmp.4, rsq = 0.6, vars = "groups")

  # Check-identicals
  # Poly-order-2
  expect_identical(test.scmp.2@Significant@genes$ColdvsControl,
    expected = sigs.2$summary$ColdvsControl
  )
  # Poly-order-3
  expect_identical(test.scmp.3@Significant@genes$ColdvsControl,
    expected = sigs.3$summary$ColdvsControl
  )
  # Poly-order-4
  expect_identical(test.scmp.4@Significant@genes$ColdvsControl,
    expected = sigs.4$summary$ColdvsControl
  )
  expect_identical(test.scmp.2@Significant@genes$HeatvsControl,
    expected = sigs.2$summary$HeatvsControl[sigs.2$summary$HeatvsControl != " "]
  )
  expect_identical(test.scmp.3@Significant@genes$HeatvsControl,
    expected = sigs.3$summary$HeatvsControl[sigs.3$summary$HeatvsControl != " "]
  )
  expect_identical(test.scmp.4@Significant@genes$HeatvsControl,
    expected = sigs.4$summary$HeatvsControl[sigs.4$summary$HeatvsControl != " "]
  )
})
