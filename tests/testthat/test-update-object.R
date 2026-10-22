# The version slot and updateObject(): objects saved by an earlier version of
# scMaSigPro are brought up to date (roadmap #17d).

current_version <- function() unname(getNamespaceVersion("scMaSigPro"))

# An object saved before the version slot existed has every other slot.
without_version <- function(ob) {
    attr(ob, "version") <- NULL
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
    expect_identical(updateObject(ob), ob)
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
    for (s in setdiff(methods::slotNames(up), "version")) {
        expect_identical(methods::slot(up, s), methods::slot(ob, s), label = s)
    }
})

test_that("updateObject() updates the SingleCellExperiment slots", {
    sce <- load_dataset("splat.sim")
    ob <- without_version(splat_scmp())
    ob@Sparse <- sce
    expect_identical(updateObject(ob)@Sparse, updateObject(sce))
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

test_that("updateObject() adds the count filter parameters (decision D17)", {
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
    expect_identical(up, ob)
    # An object saved before the version slot gets both
    older <- without_version(old)
    expect_identical(updateObject(older), ob)
})

test_that("updateObject() warns about a relabelled family (ledger LK-4)", {
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
