# Structural validity of the object and of its two SingleCellExperiment slots.
# Regression tests for roadmap items S1 and S2 (invalid Dense slot, cell
# metadata misaligned with the count matrix after binning).

test_that("objects are valid S4 at every stage of the splat.sim workflow", {
    expect_valid_scmp(splat_scmp())
    expect_valid_scmp(splat_binned())
    expect_valid_scmp(splat_design(2))
    expect_valid_scmp(splat_fitted("default"))
    expect_valid_scmp(splat_filtered())
    expect_valid_scmp(splat_clustered())
})

test_that("objects are valid S4 after sc.squeeze() and sc.restruct() on multi.lin.sce", {
    expect_valid_scmp(multi_binned())
    expect_valid_scmp(multi_restruct())
})

test_that("the bundled scmp.ob is a valid object", {
    expect_valid_scmp(load_dataset("scmp.ob"))
})

test_that("cell metadata stays aligned with the counts after binning", {
    ob <- splat_binned()
    expect_sparse_aligned(ob)

    raw <- as.matrix(SingleCellExperiment::counts(load_dataset("splat.sim")))
    via_accessor <- as.matrix(SingleCellExperiment::counts(ob@Sparse))
    expect_equal(via_accessor[, colnames(raw)], raw)
    expect_equal(eSparse(ob)[, colnames(raw)], raw)

    meta <- cSparse(ob)
    expect_identical(rownames(meta), colnames(via_accessor))
})

test_that("the Dense slot has gene rows and bin columns", {
    for (ob in list(splat_binned(), multi_binned(), multi_restruct())) {
        expect_identical(dim(ob@Dense), dim(eDense(ob)))
        expect_identical(rownames(ob@Dense), rownames(eDense(ob)))
        expect_identical(colnames(ob@Dense), rownames(cDense(ob)))
        expect_identical(colnames(eDense(ob)), rownames(cDense(ob)))
    }
})

test_that("sc.squeeze() realigns cell metadata that arrives out of order", {
    # Objects binned by 0.0.4 carry cell metadata in a different order from the
    # counts; binning such an object again must repair the pairing.
    ob <- splat_scmp()
    ob@Sparse@colData <- ob@Sparse@colData[rev(seq_len(ncol(ob@Sparse))), ]
    ob <- sc.squeeze(ob)
    expect_valid_scmp(ob)
    expect_sparse_aligned(ob)
    expect_equal(eDense(ob), eDense(splat_binned()))
})

test_that("pb_counts() names the bins after the rownames of the Dense metadata", {
    ob <- splat_binned()
    bins <- cDense(ob)
    rownames(bins) <- paste0("renamed_", rownames(bins))
    cDense(ob) <- bins
    ob <- pb_counts(ob)
    expect_valid_scmp(ob)
    expect_identical(colnames(eDense(ob)), rownames(cDense(ob)))
    expect_equal(unname(eDense(ob)), unname(eDense(splat_binned())))
    ob <- sc.set.poly(ob, poly_degree = 2)
    expect_identical(rownames(predictors(ob)), colnames(eDense(ob)))
    expect_no_error(sc.p.vector(ob, verbose = FALSE))
})

test_that("pb_counts() falls back to the bin column when bins have no rownames", {
    ob <- splat_binned()
    bins <- cDense(ob)
    rownames(bins) <- NULL
    cDense(ob) <- bins
    ob <- pb_counts(ob)
    expect_valid_scmp(ob)
    expect_identical(colnames(eDense(ob)), cDense(ob)[[ob@Parameters@bin_col]])
})

test_that("pb_counts() rejects duplicated bin names", {
    ob <- splat_binned()
    bins <- SummarizedExperiment::colData(ob@Dense)
    rownames(bins)[2] <- rownames(bins)[1]
    ob@Dense@colData <- bins
    expect_error(pb_counts(ob), "must be unique")
})

test_that("sc.squeeze() removes cells that fall outside every bin, with a warning", {
    # split_bins leaves gaps between the split bins on multi.lin.sce (3 cells
    # with Sturges); 0.0.4 dropped them from the metadata only, leaving an
    # invalid object.
    expect_warning(
        ob <- sc.squeeze(multi_scmp(), split_bins = TRUE),
        "could not be assigned to a bin"
    )
    expect_valid_scmp(ob)
    expect_sparse_aligned(ob)
    expect_identical(ncol(ob@Sparse), 1974L)
})

test_that("the prototype of every slot class is a valid object", {
    classes <- c(
        "ParameterConfig", "MatrixDesign", "VariableProfiles", "Estimates",
        "Significant"
    )
    for (cls in classes) {
        expect_true(validObject(new(cls), test = TRUE), label = cls)
    }
})

test_that("validity methods return every problem instead of stopping", {
    params <- new("ParameterConfig")
    params@path_col <- character()
    params@p_value <- c(0.05, 0.01)
    params@offset <- NA
    params@distribution <- "negative binomial"
    expect_identical(validObject(params, test = TRUE), c(
        "'path_col' must be a single non-missing value",
        "'p_value' must be a single non-missing value",
        "'offset' must be a single non-missing value",
        "'distribution' must be a 'family' object"
    ))

    design <- splat_design(2)@Design
    design@assignment_matrix <- design@assignment_matrix[-1, , drop = FALSE]
    design@offset <- design@offset[-1]
    design@groups.vector <- design@groups.vector[-1]
    expect_identical(validObject(design, test = TRUE), c(
        "'assignment_matrix' must have the rows of 'predictor_matrix'",
        "'offset' must have one value per row of 'predictor_matrix'",
        "'groups.vector' must have one value per column of 'predictor_matrix'"
    ))

    profile <- splat_fitted("default")@Profile
    profile@adj_p_values <- profile@adj_p_values[-1]
    profile@fdr <- c(0.01, 0.02)
    profile@non_flat <- c(profile@non_flat, "NotAGene")
    expect_identical(validObject(profile, test = TRUE), c(
        "'adj_p_values' must have the length and names of 'p_values'",
        "'fdr' must hold at most one value",
        "every gene in 'non_flat' must be named in 'p_values'"
    ))

    estimate <- splat_fitted("default")@Estimate
    estimate@t_score_matrix <- estimate@t_score_matrix[-1, , drop = FALSE]
    estimate@path_coefficient_matrix <-
        estimate@path_coefficient_matrix[-1, , drop = FALSE]
    rownames(estimate@significance_matrix)[1] <- "NotAGene"
    expect_identical(validObject(estimate, test = TRUE), c(
        "'t_score_matrix' must have the shape and rows of 'coefficient_matrix'",
        "'path_coefficient_matrix' must have the rows of 'coefficient_matrix'",
        "rows of 'significance_matrix' must be rows of 'coefficient_matrix'"
    ))

    significant <- splat_clustered()@Significant
    significant@genes[[1]] <- seq_len(3)
    names(significant@clusters) <- NULL
    expect_identical(validObject(significant, test = TRUE), c(
        "every element of 'genes' must be a character vector of gene names",
        "'clusters' must hold one value per gene, named by gene"
    ))
})

test_that("a new ScMaSigPro object is valid", {
    expect_true(validObject(new("ScMaSigPro"), test = TRUE))
})

test_that("validObject() names the slot of each problem in an ScMaSigPro", {
    ob <- splat_clustered()
    ob@Parameters@p_value <- c(0.05, 0.01)
    ob@Dense@colData <- ob@Dense@colData[-1, ]
    problems <- validObject(ob, test = TRUE)
    expect_length(problems, 2)
    expect_match(problems[1], "^slot 'Dense': ")
    expect_identical(
        problems[2],
        "slot 'Parameters': 'p_value' must be a single non-missing value"
    )
    expect_error(validObject(ob), "slot 'Parameters'", fixed = TRUE)
})

test_that("the steps refuse to return an invalid object", {
    # Each step ends with validObject(), directly or through pb_counts(), so
    # a problem in any slot is reported where the object is returned.
    broken <- function(ob) {
        ob@Parameters@anno_col <- character()
        ob
    }
    msg <- "slot 'Parameters': 'anno_col' must be a single non-missing value"
    expect_error(sc.squeeze(broken(splat_scmp())), msg, fixed = TRUE)
    expect_error(pb_counts(broken(splat_binned())), msg, fixed = TRUE)
    expect_error(
        sc.restruct(broken(multi_binned()),
            end_node_list = list("Y_15", "Y_18", "Y_52"),
            root_node = "Root", link_node_list = list("Y_18_links_Y_15"),
            verbose = FALSE
        ),
        msg,
        fixed = TRUE
    )
    expect_error(sc.set.poly(broken(splat_binned())), msg, fixed = TRUE)
    ob <- broken(splat_design(2))
    expect_error(sc.p.vector(ob, verbose = FALSE), msg, fixed = TRUE)
    # No gene selected: the early return is checked too
    expect_error(
        suppressMessages(sc.p.vector(ob, p_value = 0, verbose = FALSE)),
        msg,
        fixed = TRUE
    )
    ob <- broken(splat_fitted("default"))
    expect_error(sc.t.fit(ob, verbose = FALSE), msg, fixed = TRUE)
    expect_error(sc.filter(ob, vars = "groups"), msg, fixed = TRUE)
    expect_error(
        sc.cluster.trend(broken(splat_filtered()), geneSet = "union", k = 4),
        msg,
        fixed = TRUE
    )
})

test_that("arguments that would be stored as invalid parameters are errors", {
    expect_error(
        as_scmp(load_dataset("splat.sim"),
            root_label = NA_character_, verbose = FALSE,
            additional_params = list(
                labels_exist = TRUE, exist_ptime_col = "Step",
                exist_path_col = "Group"
            )
        ),
        "'root_label' must be a single character string.",
        fixed = TRUE
    )
    expect_error(
        sc.p.vector(splat_design(2), p_value = c(0.05, 0.1), verbose = FALSE),
        "slot 'Parameters': 'p_value' must be a single non-missing value",
        fixed = TRUE
    )
    expect_error(
        sc.p.vector(splat_design(2),
            offset = FALSE, log_offset = NA, verbose = FALSE
        ),
        "slot 'Parameters': 'log_offset' must be a single non-missing value",
        fixed = TRUE
    )
})

test_that("bins without counts pass the check of sc.set.poly()", {
    # cDense<- with another number of bins leaves the Dense slot without
    # counts until pb_counts(); that object is valid, so sc.set.poly() takes
    # it.
    ob <- splat_binned()
    suppressMessages(cDense(ob) <- cDense(ob)[-1, ])
    ob <- sc.set.poly(ob, poly_degree = 2)
    expect_valid_scmp(ob)
    expect_identical(nrow(predictors(ob)), 15L)
})

test_that("the count filter parameters are at least 0 (decision D17)", {
    params <- new("ParameterConfig")
    expect_identical(params@min_counts, 0)
    expect_identical(params@min_bins_expressed, 0)
    params@min_counts <- -1
    params@min_bins_expressed <- c(1, 2)
    expect_identical(validObject(params, test = TRUE), c(
        "'min_bins_expressed' must be a single non-missing value",
        "'min_counts' must be at least 0"
    ))
    params@min_bins_expressed <- -0.5
    expect_identical(validObject(params, test = TRUE), c(
        "'min_counts' must be at least 0",
        "'min_bins_expressed' must be at least 0"
    ))
})
