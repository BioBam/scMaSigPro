suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Reproduce Results of MaSigPro", {
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
    design_1 <- make.design.matrix(edesign = edesign.abiotic, degree = 1)
    design_3 <- make.design.matrix(edesign = edesign.abiotic, degree = 3)
    
    # Step-5: Set polynomial and make design
    test.scmp.1 <- sc.set.poly(test.scmp, poly_degree = 1)
    test.scmp.3 <- sc.set.poly(test.scmp, poly_degree = 3)
    
    # Check-dis
    # Order-1
    expect_equal(
        expected = design_1$dis,
        object = test.scmp.1@Design@predictor_matrix
    )
    # Poly-order-3
    expect_equal(
        expected = design_3$dis,
        object = test.scmp.3@Design@predictor_matrix
    )
    
    # Check-edesign
    # Poly-order-1
    expect_equal(
        expected = design_1$edesign,
        object = test.scmp.1@Design@assignment_matrix
    )
    # Poly-order-3
    expect_equal(
        expected = design_3$edesign,
        object = test.scmp.3@Design@assignment_matrix
    )
    # Check group-vector
    expect_identical(
        expected = design_1$groups.vector,
        object = test.scmp.1@Design@groups.vector
    )
    # Poly-order-3
    expect_identical(
        expected = design_3$groups.vector,
        object = test.scmp.3@Design@groups.vector
    )
    
    # Step-6: Run P.vector
    gc <- capture_output(fit_1 <- p.vector(data.abiotic, design_1,
                                           Q = 0.05,
                                           MT.adjust = "BH", min.obs = 20
    ))
    gc <- capture_output(fit_3 <- p.vector(data.abiotic, design_3,
                                           Q = 0.05,
                                           MT.adjust = "BH", min.obs = 20
    ))
    
    # Step-7: Run sc.p.vector
    test.scmp.1 <- sc.p.vector(test.scmp.1,
                               min_na = 20, verbose = FALSE,
                               offset = FALSE, parallel = FALSE, max_it = 25,
                               epsilon = 0.00001, family = gaussian()
    )
    test.scmp.3 <- sc.p.vector(test.scmp.3,
                               min_na = 20, verbose = FALSE,
                               offset = FALSE, parallel = FALSE, max_it = 25,
                               epsilon = 0.00001, family = gaussian()
    )
    
    # Check-fdr
    # Poly-order-1
    expect_identical(
        expected = fit_1$FDR,
        object = as.vector(test.scmp.1@Profile@fdr)
    )
    # Poly-order-3
    expect_identical(
        expected = fit_3$FDR,
        object = as.vector(test.scmp.3@Profile@fdr)
    )
    
    # Check adjusted p-values
    # Poly-order-2
    expect_identical(
        expected = as.vector(fit_1$p.adjusted),
        object = as.vector(test.scmp.1@Profile@adj_p_values)
    )
    # Poly-order-3
    expect_identical(
        expected = as.vector(fit_3$p.adjusted),
        object = as.vector(test.scmp.3@Profile@adj_p_values)
    )
    
    # Check-p.value
    # Poly-order-1
    expect_identical(
        expected = as.vector(fit_1$p.vector),
        object = as.vector(test.scmp.1@Profile@p_values)
    )
    # Poly-order-3
    expect_identical(
        expected = as.vector(fit_3$p.vector),
        object = as.vector(test.scmp.3@Profile@p_values)
    )
    
    # Check non-flat genes
    # Poly-order-2
    expect_identical(
        expected = rownames(fit_1$SELEC),
        object = test.scmp.1@Profile@non_flat
    )
    # Poly-order-3
    expect_identical(
        expected = rownames(fit_3$SELEC),
        object = test.scmp.3@Profile@non_flat
    )
    
    # Step-8: Run t.fit
    gc <- capture_output(tstep_1 <- T.fit(fit_1, step.method = "backward", alfa = 0.05))
    gc <- capture_output(tstep_3 <- T.fit(fit_3, step.method = "backward", alfa = 0.05))
    
    # Step-9: Run sc.t.fit
    test.scmp.1 <- sc.t.fit(test.scmp.1, verbose = FALSE)
    test.scmp.3 <- sc.t.fit(test.scmp.3, verbose = FALSE)
    
    # Check-tscore
    # Poly-order-1
    expect_identical(as.data.frame(showTS(test.scmp.1)), expected = tstep_1$t.score)
    # Poly-order-3
    expect_identical(as.data.frame(showTS(test.scmp.3)), expected = tstep_3$t.score)
    
    # Check-sol
    # Poly-order-1
    expect_identical(showSol(test.scmp.1), expected = tstep_1$sol)
    # Poly-order-3
    expect_identical(showSol(test.scmp.3), expected = tstep_3$sol)
    
    # Check group-coefficients
    # Poly-order-1
    expect_identical(as.data.frame(showGroupCoeff(test.scmp.1)), expected = tstep_1$group.coeffs)
    # Poly-order-3
    expect_identical(as.data.frame(showGroupCoeff(test.scmp.3)), expected = tstep_3$group.coeffs)
    
    # Check coefficients
    # Poly-order-1
    expect_identical(showCoeff(test.scmp.1), expected = tstep_1$coefficients)
    # Poly-order-3
    expect_identical(showCoeff(test.scmp.3), expected = tstep_3$coefficients)
    
    # Check influ
    # Poly-order-1
    expect_identical(showInflu(test.scmp.1), expected = tstep_1$influ.info)
    # Poly-order-3
    expect_identical(showInflu(test.scmp.3), expected = tstep_3$influ.info)
    
    # Step-10: Compute Clusters with maSigPro
    sigs.1 <- get.siggenes(tstep_1, rsq = 0.6, vars = "groups")
    sigs.3 <- get.siggenes(tstep_3, rsq = 0.6, vars = "groups")
    
    # Step-11: Compute Clusters with scMaSigPro
    test.scmp.1 <- sc.filter(test.scmp.1, rsq = 0.6, vars = "groups")
    test.scmp.3 <- sc.filter(test.scmp.3, rsq = 0.6, vars = "groups")
    
    # Check-identicals
    # Poly-order-2
    expect_identical(test.scmp.1@Significant@genes$ColdvsControl,
                     expected = sigs.1$summary$ColdvsControl
    )
    # Poly-order-3
    expect_identical(test.scmp.3@Significant@genes$ColdvsControl,
                     expected = sigs.3$summary$ColdvsControl
    )
    
    expect_identical(test.scmp.1@Significant@genes$HeatvsControl,
                     expected = sigs.1$summary$HeatvsControl[sigs.1$summary$HeatvsControl != " "]
    )
    expect_identical(test.scmp.3@Significant@genes$HeatvsControl,
                     expected = sigs.3$summary$HeatvsControl[sigs.3$summary$HeatvsControl != " "]
    )
})
