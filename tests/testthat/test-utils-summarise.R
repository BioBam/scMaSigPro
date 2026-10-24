# .summarise_by() replaces dplyr::group_by() + summarise(.groups = "drop")
# in plotTrend() and plotTrendCluster().

test_that(".summarise_by() sorts groups as dplyr did and keeps row order", {
    df <- data.frame(
        g = c("b", "B", "a", "b", NA, "a"),
        x = c(2, 1, 2, 2, 1, 1),
        v = c(10, 20, 30, 40, 50, 60)
    )
    first <- .summarise_by(df, c("x", "g"), "v", function(v) v[1])
    expect_identical(colnames(first), c("x", "g", "v"))
    expect_identical(first$x, c(1, 1, 1, 2, 2))
    # C-locale order (upper case first), NA as the last group
    expect_identical(first$g, c("B", "a", NA, "a", "b"))
    # the first value of each group, in row order
    expect_identical(first$v, c(20, 60, 50, 30, 10))
    expect_identical(.row_names_info(first), -5L)
    expect_identical(.summarise_by(df, "g", "v", sum)$v, c(20, 90, 50, 50))
})
