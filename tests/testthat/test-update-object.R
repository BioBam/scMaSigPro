# The version slot and updateObject(): objects saved by an earlier version of
# scMaSigPro are brought up to date.

current_version <- function() unname(getNamespaceVersion("scMaSigPro"))

# An object saved before the version slot existed has every other slot.
without_version <- function(ob) {
    attr(ob, "version") <- NULL
    ob
}

# The object as updateObject() returns it when it has nothing of its own to
# update: the Sparse and Dense slots passed through the updateObject() of
# SingleCellExperiment, which stamps them with the installed version of
# SingleCellExperiment when they were saved by an older one (as the bundled
# scmp.ob is under a newer SingleCellExperiment)
with_current_sce <- function(ob) {
    ob@Sparse <- updateObject(ob@Sparse)
    ob@Dense <- updateObject(ob@Dense)
    ob
}

test_that("new objects carry the version of the package", {
    expect_identical(new("ScMaSigPro")@version, current_version())
    expect_identical(splat_scmp()@version, current_version())
    expect_identical(multi_scmp()@version, current_version())
    # The steps keep it
    expect_identical(splat_clustered()@version, current_version())
    expect_identical(multi_restruct()@version, current_version())
})

test_that("the version slot holds a single version number", {
    ob <- splat_scmp()
    for (v in list(character(), c("0.0.5", "0.0.6"), NA_character_, "x")) {
        ob@version <- v
        expect_identical(
            validObject(ob, test = TRUE),
            "'version' must be a single version number"
        )
    }
})

test_that("updateObject() returns a current object unchanged", {
    path <- tempfile(fileext = ".rds")
    saveRDS(load_dataset("scmp.ob"), path)
    ob <- readRDS(path)
    unlink(path)
    expect_true(validObject(ob, test = TRUE))
    expect_identical(updateObject(ob), with_current_sce(ob))
    # Objects made by the installed packages carry their versions
    expect_identical(updateObject(splat_clustered()), splat_clustered())
    expect_identical(updateObject(multi_restruct()), multi_restruct())
})

test_that("updateObject() adds the version slot to an object saved without", {
    ob <- load_dataset("scmp.ob")
    old <- without_version(ob)
    expect_error(validObject(old), "but not in object: \"version\"",
        fixed = TRUE
    )
    expect_output(print(old), "update it with updateObject()", fixed = TRUE)
    msgs <- capture_messages(up <- updateObject(old, verbose = TRUE))
    expect_match(msgs, "Adding the 'version' slot", fixed = TRUE, all = FALSE)
    expect_valid_scmp(up)
    expect_identical(up@version, current_version())
    expected <- with_current_sce(ob)
    for (s in setdiff(methods::slotNames(up), "version")) {
        expect_identical(methods::slot(up, s), methods::slot(expected, s),
            label = s
        )
    }
})

test_that("updateObject() updates the SingleCellExperiment slots", {
    sce <- load_dataset("splat.sim")
    ob <- without_version(splat_scmp())
    ob@Sparse <- sce
    expect_identical(updateObject(ob)@Sparse, updateObject(sce))
})

test_that("updateObject() updates the stamp of older SingleCellExperiments", {
    # As the bundled scmp.ob looks to a newer SingleCellExperiment
    ob <- load_dataset("scmp.ob")
    stale <- ob
    old_version <- package_version("1.11.3")
    SingleCellExperiment::int_metadata(stale@Sparse)$version <- old_version
    SingleCellExperiment::int_metadata(stale@Dense)$version <- old_version
    up <- updateObject(stale)
    stamp <- function(sce) SingleCellExperiment::int_metadata(sce)$version
    installed <- utils::packageVersion("SingleCellExperiment")
    expect_identical(stamp(up@Sparse), installed)
    expect_identical(stamp(up@Dense), installed)
    expect_identical(up, with_current_sce(stale))
})

test_that("updateObject() refuses an object binned by 0.0.4", {
    # 0.0.4 saved a Dense slot that is not a valid SingleCellExperiment.
    old_ob <- without_version(load_dataset("scmp.ob"))
    old_ob@Dense@colData <- old_ob@Dense@colData[-1, ]
    expect_error(updateObject(old_ob), "binned by scMaSigPro 0.0.4",
        fixed = TRUE
    )
    cmd <- paste0(
        "new_ob <- sc.squeeze(as_scmp(old_ob@Sparse, ",
        "ptime_col = \"Pseudotime\", path_col = \"Path\"))"
    )
    expect_error(updateObject(old_ob), cmd, fixed = TRUE)
    # The command in the message bins the cells again
    new_ob <- suppressMessages(eval(parse(text = cmd)))
    expect_valid_scmp(new_ob)
    expect_identical(eDense(new_ob), eDense(splat_binned()))
})

test_that("as_scmp() stores an up-to-date SingleCellExperiment", {
    # splat.sim was saved by an older SingleCellExperiment. colData<- updates
    # it when labels_exist = TRUE; as_scmp() now does so on both routes.
    sce <- load_dataset("splat.sim")
    ob <- as_scmp(sce, verbose = FALSE)
    expect_identical(ob@Sparse, updateObject(sce))
    expect_identical(updateObject(ob), ob)
})

test_that("updateObject() adds the count filter parameters", {
    ob <- load_dataset("scmp.ob")
    old <- ob
    params <- old@Parameters
    attr(params, "min_counts") <- NULL
    attr(params, "min_bins_expressed") <- NULL
    old@Parameters <- params
    old@version <- "0.6.0"
    expect_error(validObject(old), "but not in object", fixed = TRUE)
    msgs <- capture_messages(up <- updateObject(old, verbose = TRUE))
    expect_match(msgs,
        "Adding the parameters min_counts and min_bins_expressed (0).",
        fixed = TRUE, all = FALSE
    )
    expect_valid_scmp(up)
    expect_identical(up@Parameters@min_counts, 0)
    expect_identical(up@Parameters@min_bins_expressed, 0)
    expect_identical(up@version, current_version())
    expect_identical(up, with_current_sce(ob))
    # An object saved before the version slot gets both
    older <- without_version(old)
    expect_identical(updateObject(older), with_current_sce(ob))
})

test_that("updateObject() adds the test parameter", {
    ob <- load_dataset("scmp.ob")
    old <- ob
    params <- old@Parameters
    attr(params, "test") <- NULL
    old@Parameters <- params
    expect_error(validObject(old), "but not in object", fixed = TRUE)
    msgs <- capture_messages(up <- updateObject(old, verbose = TRUE))
    expect_match(msgs, "Adding the parameter test (\"chisq\").",
        fixed = TRUE, all = FALSE
    )
    expect_valid_scmp(up)
    expect_identical(up@Parameters@test, "chisq")
    expect_identical(up, with_current_sce(ob))
    # An object saved before the count filter and the version slot
    attr(params, "min_counts") <- NULL
    attr(params, "min_bins_expressed") <- NULL
    old@Parameters <- params
    expect_identical(updateObject(without_version(old)), with_current_sce(ob))
})

test_that("updateObject() warns about a relabelled family", {
    old <- splat_relabelled()
    msg <- paste(
        "'family' is labelled with the 'identity' link but has the link",
        "functions of the 'log' link."
    )
    expect_warning(up <- updateObject(old), msg, fixed = TRUE)
    expect_identical(up, old)
    expect_warning(updateObject(without_version(old)), "run sc.p.vector()",
        fixed = TRUE
    )
    # sc.p.vector() again gives a family with the identity link
    refit <- suppressWarnings(sc.p.vector(old,
        link = "identity", offset = FALSE, verbose = FALSE
    ))
    expect_no_warning(updateObject(refit))
    expect_no_warning(updateObject(splat_clustered()))
})

# An object as an earlier version saved it after sc.p.vector() and
# sc.t.fit(), with the log_offset parameter. With log_offset = FALSE the
# offsets of its fits were the size factors themselves
with_log_offset <- function(ob, log_offset) {
    params <- ob@Parameters
    attr(params, "log_offset") <- log_offset
    ob@Parameters <- params
    if (!log_offset && ob@Parameters@offset) {
        ob@Design@offset <- exp(ob@Design@offset)
    }
    ob
}

test_that("updateObject() removes the log_offset parameter", {
    ob <- splat_fitted("default")
    old <- with_log_offset(ob, TRUE)
    expect_true(validObject(old, test = TRUE))
    expect_error(sc.t.fit(old, verbose = FALSE),
        "update it with updateObject() first",
        fixed = TRUE
    )
    expect_output(print(old), "update it with updateObject()", fixed = TRUE)
    msgs <- capture_messages(up <- updateObject(old, verbose = TRUE))
    expect_match(msgs, "Removing the parameter log_offset.",
        fixed = TRUE, all = FALSE
    )
    # Offsets saved on the log scale are kept
    expect_identical(up, ob)
    # A fit without offset: nothing else to update
    no_offset <- splat_fitted("no_offset")
    expect_no_warning(up <- updateObject(with_log_offset(no_offset, FALSE)))
    expect_identical(up, no_offset)
})

test_that("updateObject() stores the log of size factors saved as offsets", {
    ob <- splat_fitted("default")
    old <- with_log_offset(ob, FALSE)
    expect_warning(up <- updateObject(old),
        "run sc.p.vector() and the later steps on the object again.",
        fixed = TRUE
    )
    expect_valid_scmp(up)
    expect_false(methods::.hasSlot(up@Parameters, "log_offset"))
    expect_identical(up@Design@offset, log(old@Design@offset))
    expect_equal(up@Design@offset, ob@Design@offset, tolerance = 1e-12)
    # The fits are kept as they are
    expect_identical(up@Profile, old@Profile)
    expect_identical(up@Estimate, old@Estimate)
    # An object saved before the version slot is updated the same way
    expect_warning(older <- updateObject(without_version(old)),
        "used the size factors of the bins as the offset",
        fixed = TRUE
    )
    expect_identical(older, up)
})

test_that("updateObject() removes the offsets of fits without offset", {
    # Earlier versions stored zeros in sc.set.poly() and kept the offsets of
    # an earlier fit after a fit without offset
    design <- splat_design(2)
    old <- with_log_offset(design, FALSE)
    old@Design@offset <- stats::setNames(
        numeric(nrow(predictors(design))), rownames(predictors(design))
    )
    expect_no_warning(up <- updateObject(old))
    expect_identical(up, design)
    fitted <- splat_fitted("no_offset")
    stale <- with_log_offset(fitted, FALSE)
    stale@Design@offset <- exp(splat_fitted("default")@Design@offset)
    expect_no_warning(up <- updateObject(stale))
    expect_identical(up, fitted)
})

test_that("updateObject() warns about results computed from other slots", {
    # Objects as earlier versions saved them after a step was run again:
    # the slots of the later steps kept the results computed from the
    # slots that the step replaced
    ob <- splat_clustered()
    sol <- showSol(ob)
    lists <- unique(unlist(ob@Significant@genes))
    stale <- list()
    stale$rebinned <- ob
    stale$rebinned@Dense <- sc.squeeze(splat_scmp(), bin_method = "Sqrt")@Dense
    stale$redesigned <- ob
    stale$redesigned@Design <- splat_design(3)@Design
    stale$new_terms <- stale$redesigned
    stale$new_terms@Design@offset <- ob@Design@offset
    stale$retested <- ob
    stale$retested@Profile <- sc.p.vector(splat_design(2),
        p_value = 0.001, verbose = FALSE
    )@Profile
    stale$relisted <- ob
    stale$relisted@Significant@genes$extra <- setdiff(
        rownames(ob@Sparse), rownames(sol)
    )[1]
    stale$reclustered <- ob
    stale$reclustered@Significant@clusters <- c(
        ob@Significant@clusters,
        stats::setNames(list(1L), setdiff(rownames(sol), lists)[1])
    )
    reasons <- c(
        rebinned = "its design was built for other bins",
        redesigned = "its tests used an offset that its design does not hold",
        new_terms = "its models are those of other genes or of another design",
        retested = "its models are those of other genes or of another design",
        relisted = "its gene lists hold genes without a model",
        reclustered = "its clusters hold genes that are not in its gene lists"
    )
    steps <- c(
        rebinned = "sc.set.poly()", redesigned = "sc.p.vector()",
        new_terms = "sc.t.fit()", retested = "sc.t.fit()",
        relisted = "sc.filter()", reclustered = "sc.cluster.trend()"
    )
    for (n in names(stale)) {
        expect_warning(up <- updateObject(stale[[n]]),
            paste0(
                "Some results of this object were not computed from its ",
                "earlier slots: ", reasons[[n]], ". Earlier versions of ",
                "scMaSigPro kept the results of the later steps when a step ",
                "was run again: run ", steps[[n]], " and the later steps on ",
                "the object again."
            ),
            fixed = TRUE, info = n
        )
        # The object is returned as it is
        expect_identical(up, stale[[n]], info = n)
    }
    # Controls: objects of the workflow, from every step, and bins renamed
    # after the design, which the fits accept
    renamed <- splat_design(2)
    bins <- cDense(renamed)
    rownames(bins) <- paste0("bin", seq_len(nrow(bins)))
    cDense(renamed) <- bins
    controls <- list(
        load_dataset("scmp.ob"), splat_scmp(), splat_binned(), splat_design(2),
        splat_fitted("default"), splat_fitted("no_offset"), splat_filtered(),
        splat_clustered(), multi_restruct(), renamed
    )
    for (i in seq_along(controls)) {
        expect_no_warning(updateObject(controls[[i]]))
    }
})
