# Functions that use a suggested package check that it is installed and call
# it with `::`.

test_that("plotIntersect() draws an UpSetR plot", {
    skip_if_not_installed("UpSetR")
    # UpSetR uses ggplot2 features deprecated upstream; silence their
    # lifecycle warnings, which are not about scMaSigPro
    op <- options(lifecycle_verbosity = "quiet")
    on.exit(options(op), add = TRUE)
    # upset() builds grobs, which opens a graphics device: use one that
    # writes no Rplots.pdf
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_s3_class(plotIntersect(splat_filtered()), "upset")
})

test_that("Mclust clustering needs mclust attached (ledger M1)", {
    skip_if_not_installed("mclust")
    ob <- splat_filtered()
    run <- function() {
        sc.cluster.trend(ob,
            geneSet = "union", k = 4, cluster_method = "Mclust"
        )
    }
    # Mclust() evaluates mclustBIC() in its caller's frame, where it is not
    # visible unless mclust is attached; frozen as it is until it is fixed
    if (!"package:mclust" %in% search()) {
        expect_error(
            run(), "could not find function \"mclustBIC\"", fixed = TRUE
        )
        suppressPackageStartupMessages(attachNamespace("mclust"))
        on.exit(detach("package:mclust"), add = TRUE)
    }
    clusters <- run()@Significant@clusters
    expect_setequal(names(clusters), Reduce(union, ob@Significant@genes))
    expect_true(all(unlist(clusters) %in% 1:4))
})
