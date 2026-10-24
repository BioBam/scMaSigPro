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

test_that("Mclust clustering works without attaching mclust (ledger M1)", {
    skip_if_not_installed("mclust")
    ob <- splat_filtered()
    run <- function() {
        sc.cluster.trend(ob,
            geneSet = "union", k = 4, cluster_method = "Mclust"
        )@Significant@clusters
    }
    # Mclust() evaluates mclustBIC() in the frame it is called from; the
    # clusters must not depend on whether mclust is attached
    if (!"package:mclust" %in% search()) {
        clusters <- run()
        suppressPackageStartupMessages(attachNamespace("mclust"))
        on.exit(detach("package:mclust"), add = TRUE)
        expect_identical(run(), clusters)
    }
    clusters <- run()
    expect_setequal(names(clusters), Reduce(union, ob@Significant@genes))
    expect_true(all(unlist(clusters) %in% 1:4))
})
