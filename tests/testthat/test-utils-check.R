# .check() replaces assertthat::assert_that(cond, msg = ...): it passes only
# a single TRUE and stops with the message alone (no call).

test_that(".check() passes a single TRUE invisibly", {
    expect_invisible(.check(TRUE, "unused"))
    expect_true(.check(c(a = TRUE), "unused"))
})

test_that(".check() stops with its message for anything else", {
    for (cond in list(FALSE, NA, logical(0), c(TRUE, TRUE), 1, "TRUE", NULL)) {
        expect_error(.check(cond, "custom message"), "custom message",
            fixed = TRUE
        )
    }
})

test_that(".check() errors without a call and evaluates msg lazily", {
    err <- tryCatch(.check(FALSE, "no call"), error = identity)
    expect_null(conditionCall(err))
    expect_identical(conditionMessage(err), "no call")
    err <- tryCatch(.check(FALSE, c("first.", "second.")), error = identity)
    expect_identical(conditionMessage(err), "first.\nsecond.")
    expect_silent(.check(TRUE, stop("msg must not be evaluated")))
})
