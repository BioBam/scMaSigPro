# Every gene is fitted on its own with the same glm() calls whatever BPPARAM
# runs them, so a parallel backend must reproduce the serial fit exactly.

expect_serial_fit <- function(bpparam) {
    ref <- splat_fitted("default")
    pv <- sc.p.vector(splat_design(2), verbose = FALSE, BPPARAM = bpparam)
    expect_identical(pv@Profile, ref@Profile)
    tf <- sc.t.fit(pv, verbose = FALSE, BPPARAM = bpparam)
    expect_identical(tf@Estimate, ref@Estimate)
    expect_identical(tf@Design, ref@Design)
}

# SnowParam workers are new R sessions that look the package namespace up by
# name when they unserialise the per-gene function. Under R CMD check they find
# the installed package being checked; under load_all() (devtools::test()) they
# would find an installed copy, if any, instead of the source under test, so
# each worker loads that source first. Returns the namespace path of each
# distinct worker.
snow_namespace_paths <- function(p) {
    on_worker <- function(i, path, dev) {
        if (dev) {
            pkgload::load_all(path,
                export_all = FALSE, helpers = FALSE,
                attach_testthat = FALSE, quiet = TRUE
            )
        }
        ns_path <- getNamespaceInfo("scMaSigPro", "path")
        c(pid = Sys.getpid(), path = normalizePath(ns_path))
    }
    # A function of the global environment does not make the worker load the
    # namespace before load_all() runs.
    environment(on_worker) <- globalenv()
    res <- BiocParallel::bplapply(
        seq_len(BiocParallel::bpnworkers(p)), on_worker,
        path = getNamespaceInfo("scMaSigPro", "path"),
        dev = pkgload::is_dev_package("scMaSigPro"), BPPARAM = p
    )
    res <- do.call(rbind, res)
    unname(res[!duplicated(res[, "pid"]), "path"])
}

# BiocParallel picks a socket port from the clock and an RNG stream that
# starts from the same seed in every R session, so two test runs on one
# machine can pick the same port. Derive the port from the process id instead,
# and try the next one if it cannot be opened.
start_param <- function(make_param) {
    for (port in 11000L + (Sys.getpid() + 0:4) %% 1000L) {
        p <- tryCatch(BiocParallel::bpstart(make_param(port)),
            warning = function(w) NULL, error = function(e) NULL
        )
        if (!is.null(p)) {
            return(p)
        }
    }
    stop("could not start a BiocParallel cluster")
}

test_that("SerialParam() is the default and is reported when verbose", {
    expect_identical(formals(sc.p.vector)$BPPARAM, quote(SerialParam()))
    expect_identical(formals(sc.t.fit)$BPPARAM, quote(SerialParam()))
    expect_message(
        sc.p.vector(splat_design(2)),
        "Fitting the full model to 100 genes with SerialParam (workers: 1)",
        fixed = TRUE
    )
    expect_serial_fit(BiocParallel::SerialParam())
})

test_that("MulticoreParam(2) reproduces the serial fit", {
    skip_on_os("windows")
    p <- start_param(function(port) {
        BiocParallel::MulticoreParam(2, manager.port = port)
    })
    on.exit(BiocParallel::bpstop(p), add = TRUE)
    expect_serial_fit(p)
})

test_that("SnowParam(2) reproduces the serial fit", {
    skip_on_cran()
    # exportglobals = FALSE: under testthat, options() hold the test
    # environment, and exporting them would make every worker load an
    # installed copy of the package while it unserialises them, before
    # snow_namespace_paths() runs.
    p <- start_param(function(port) {
        BiocParallel::SnowParam(2, exportglobals = FALSE, manager.port = port)
    })
    on.exit(BiocParallel::bpstop(p), add = TRUE)
    path <- normalizePath(getNamespaceInfo("scMaSigPro", "path"))
    expect_identical(snow_namespace_paths(p), c(path, path))
    expect_serial_fit(p)
})
