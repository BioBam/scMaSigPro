suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(scMaSigPro))
suppressPackageStartupMessages(library(maSigPro))

test_that("Check-'design$groups.vector' Reproducibility: Match, Dimension, Name and Value", {
    
    # Step-1: Load Data
    data("data.abiotic")
    data("edesign.abiotic")
    
    # Create MaSigPro Design
    design_2 <- make.design.matrix(edesign = edesign.abiotic, degree = 2)
    design_3 <- make.design.matrix(edesign = edesign.abiotic, degree = 3)
    design_4 <- make.design.matrix(edesign = edesign.abiotic, degree = 4)
    
    # Convert the Design for scMaSigPro
    add_col <- function(x){
        return(names(x[x==1]))
    }
    edesign.abiotic.in <- as.data.frame(edesign.abiotic[, 1, drop = F])
    Group <- c(unlist(apply(edesign.abiotic[, c(3:6)], 1, add_col, simplify = F), use.names = F))
    edesign.abiotic.in$Group <- Group
    
    # Create Design with scMaSigPro
    test.scmp <- create_scmpObj(counts = data.abiotic,
                                cell_data = edesign.abiotic.in,
                                pseudotime_colname = "Time",
                                path_colname = "Group",
                                use_as_bin = T)
    
    # Run Make Design
    test.scmp.2 <- sc.make.design.matrix(test.scmp, poly_degree = 2)
    test.scmp.3 <- sc.make.design.matrix(test.scmp, poly_degree = 3)
    test.scmp.4 <- sc.make.design.matrix(test.scmp, poly_degree = 4)
    
    # Check
    # Poly-order-2
    expect_equal(expected = design_2$groups.vector,
                 object = test.scmp.2@edesign@groups.vector)
    
    # Poly-order-3
    expect_equal(expected = design_3$groups.vector,
                 object = test.scmp.3@edesign@groups.vector)
    
    # Poly-order-4
    expect_equal(expected = design_4$groups.vector,
                 object = test.scmp.4@edesign@groups.vector)
})