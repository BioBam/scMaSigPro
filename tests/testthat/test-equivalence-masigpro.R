# Oracle tests: with no offset, scMaSigPro must reproduce maSigPro exactly.
# test-reproduce-masigpro.R covers the Gaussian/identity case on the full
# data.abiotic with backward selection. This file extends the oracle to every
# selection method and polynomial degree (on a 200-gene subset, for speed), and
# to the negative binomial family on maSigPro's NBdata.

suppressPackageStartupMessages(library(maSigPro))

quiet <- function(expr) {
  invisible(utils::capture.output(res <- expr))
  res
}

group_column <- function(design, group_cols) {
  md <- as.data.frame(design)
  md$Group <- apply(md[, group_cols, drop = FALSE], 1, function(x) names(x)[x == 1])
  md[, !(colnames(md) %in% group_cols), drop = FALSE]
}

methods <- c("backward", "forward", "two.ways.backward", "two.ways.forward")

# ---- Gaussian, data.abiotic subset -------------------------------------------

data("data.abiotic", package = "maSigPro")
data("edesign.abiotic", package = "maSigPro")
abiotic <- as.matrix(data.abiotic[seq_len(200), ])
abiotic_meta <- group_column(edesign.abiotic, c("Control", "Cold", "Heat", "Salt"))

abiotic_scmp <- create_scmp(
  counts = abiotic, cell_data = abiotic_meta,
  ptime_col = "Time", path_col = "Group", use_as_bin = TRUE
)

for (degree in 1:3) {
  design <- make.design.matrix(edesign.abiotic, degree = degree)
  fit <- quiet(p.vector(abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20))

  ob <- sc.set.poly(abiotic_scmp, poly_degree = degree)
  ob <- sc.p.vector(ob,
    min_na = 20, verbose = FALSE, link = "identity", offset = FALSE,
    parallel = FALSE, max_it = 25, epsilon = 0.00001, family = gaussian()
  )

  test_that(paste0("Gaussian degree ", degree, ": design and p.vector match maSigPro"), {
    expect_equal(predictors(ob), as.matrix(design$dis))
    expect_equal(pathAssign(ob), as.matrix(design$edesign))
    expect_identical(ob@Design@groups.vector, design$groups.vector)
    expect_identical(as.vector(ob@Profile@p_values), as.vector(fit$p.vector))
    expect_identical(as.vector(ob@Profile@adj_p_values), as.vector(fit$p.adjusted))
    expect_identical(ob@Profile@non_flat, rownames(fit$SELEC))
  })

  for (method in methods) {
    test_that(paste0("Gaussian degree ", degree, ", ", method, ": T.fit matches maSigPro"), {
      tstep <- quiet(T.fit(fit, step.method = method, alfa = 0.05))
      sc <- sc.t.fit(ob, selection_method = method, verbose = FALSE)
      expect_equal(showSol(sc), tstep$sol)
      expect_equal(showCoeff(sc), tstep$coefficients)
      expect_equal(as.data.frame(showGroupCoeff(sc)), tstep$group.coeffs)
      expect_equal(as.data.frame(showTS(sc)), tstep$t.score)
    })
  }
}

# ---- Negative binomial, NBdata (theta = 10, log link, no offset) ------------

data("NBdata", package = "maSigPro")
data("NBdesign", package = "maSigPro")
nb_meta <- group_column(NBdesign, c("Group.1", "Group.2"))

nb_scmp <- create_scmp(
  counts = as.matrix(NBdata), cell_data = nb_meta,
  ptime_col = "Time", path_col = "Group", use_as_bin = TRUE
)
nb_scmp <- sc.set.poly(nb_scmp, poly_degree = 2)
nb_design <- make.design.matrix(NBdesign, degree = 2)

nb_fit <- quiet(p.vector(NBdata, nb_design,
  Q = 0.05, MT.adjust = "BH", min.obs = 6, counts = TRUE, theta = 10
))
nb_scmp <- sc.p.vector(nb_scmp,
  min_na = 6, verbose = FALSE, link = "log", offset = FALSE,
  family = MASS::negative.binomial(theta = 10), epsilon = 1e-5, max_it = 25
)

test_that("negative binomial p.vector matches maSigPro (counts = TRUE, theta = 10)", {
  expect_identical(as.vector(nb_scmp@Profile@p_values), as.vector(nb_fit$p.vector))
  expect_identical(nb_scmp@Profile@non_flat, rownames(nb_fit$SELEC))
})

for (method in methods) {
  test_that(paste0("negative binomial, ", method, ": T.fit matches maSigPro"), {
    tstep <- quiet(T.fit(nb_fit,
      step.method = method, alfa = 0.05,
      family = MASS::negative.binomial(theta = 10), epsilon = 1e-5
    ))
    sc <- suppressMessages(sc.t.fit(nb_scmp, selection_method = method, verbose = FALSE))
    expect_equal(showSol(sc), tstep$sol)
    expect_equal(showCoeff(sc), tstep$coefficients)
  })
}
