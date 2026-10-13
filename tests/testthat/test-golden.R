# Golden tests: the outputs of the whole workflow on the bundled datasets,
# frozen at version 0.0.4. Any difference means a change in results; if the
# change is intended, regenerate the fixture with fixtures/make-golden.R and
# document it in NEWS.

golden <- readRDS(test_path("fixtures", "golden-v0.0.4.rds"))
current <- build_golden()

# Numerical tolerance for fitted quantities (p-values, coefficients): large
# enough to absorb BLAS differences between platforms, small enough that any
# real change in the model is caught. Discrete results compare exactly.
TOL <- 1e-6

# waldo (like all.equal) switches to an absolute difference when the mean
# reference value is below the tolerance, so on the raw scale every p-value
# below 1e-6 is effectively unchecked. Compare p-value leaves on the log10
# scale as well; 1e-300 floors exact zeros (underflow).
is_p_leaf <- function(path) grepl("^fits/[^/]+/(p_values|adj_p_values|sol)$", path)
log_p <- function(x) log10(pmax(x, 1e-300))

walk_compare <- function(cur, gold, path) {
    if (is.list(gold) && !is.data.frame(gold) && !is.null(names(gold))) {
        expect_identical(names(cur), names(gold), label = paste(path, "names"))
        for (nm in names(gold)) {
            walk_compare(cur[[nm]], gold[[nm]], paste(path, nm, sep = "/"))
        }
    } else if (is.list(gold) && !is.data.frame(gold)) {
        expect_identical(length(cur), length(gold), label = paste(path, "length"))
        for (i in seq_along(gold)) {
            walk_compare(cur[[i]], gold[[i]], paste0(path, "[[", i, "]]"))
        }
    } else {
        expect_equal(cur, gold, tolerance = TOL, label = path, expected.label = "golden")
        if (is_p_leaf(path)) {
            expect_equal(log_p(cur), log_p(gold),
                tolerance = TOL,
                label = paste(path, "(log10)"), expected.label = "golden"
            )
        }
    }
}

for (section in setdiff(names(golden), "meta")) {
    test_that(paste("golden:", section), {
        expect_true(section %in% names(current))
        walk_compare(current[[section]], golden[[section]], section)
    })
}
