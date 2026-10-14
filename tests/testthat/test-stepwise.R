# The four stepwise selectors of R/stepwise.R on small problems whose answer
# can be worked out by hand. test-equivalence-masigpro.R and
# test-reproduce-masigpro.R check them against maSigPro on real data.

selectors <- list(
    backward = sc.stepback,
    forward = sc.stepfor,
    two.ways.backward = sc.two.ways.stepback,
    two.ways.forward = sc.two.ways.stepfor
)

select_with <- function(f, y, d, ...) {
    f(
        y = y, d = d, alfa = 0.05, family = gaussian(), epsilon = 1e-5,
        useOffset = NULL, useWeight = NULL, max_it = 25, ...
    )
}

# Eight observations. x1 is a trend, x2 and x3 are contrasts orthogonal to
# the intercept, to x1 and to each other, and e is orthogonal to all of them.
# So in y = 1 + 2 * x1 + e the coefficients of x2 and x3 are exactly 0
# (p-value 1) and x1 is very significant, in every submodel.
x1 <- 1:8
x2 <- c(1, -1, 1, -1, -1, 1, -1, 1)
x3 <- c(1, 1, -1, -1, -1, -1, 1, 1)
e <- c(0.2, -0.2, -0.2, 0.2, -0.2, 0.2, 0.2, -0.2)
d <- data.frame(x1 = x1, x2 = x2, x3 = x3, row.names = paste0("bin", 1:8))

test_that("the design of the hand-made problems is orthogonal", {
    cross <- crossprod(cbind(1, x1 - mean(x1), x2, x3, e))
    expect_equal(cross[upper.tri(cross)], rep(0, 10))
})

test_that("every method keeps the significant term and drops the noise", {
    y <- 1 + 2 * x1 + e
    for (method in names(selectors)) {
        fit <- select_with(selectors[[method]], y, d)
        expect_identical(names(coef(fit)), c("(Intercept)", "x1"),
            info = method
        )
        expect_equal(unname(coef(fit)), c(1, 2), info = method)
    }
})

test_that("every method returns the intercept-only model for pure noise", {
    y <- 5 + e
    for (method in names(selectors)) {
        fit <- select_with(selectors[[method]], y, d)
        expect_identical(names(coef(fit)), "(Intercept)", info = method)
        expect_equal(unname(coef(fit)), 5, info = method)
        expect_identical(attr(terms(fit), "term.labels"), character(0))
    }
})

test_that("forward methods list the terms in the order they entered", {
    # Alone, x3 explains far more than x1, so it enters first
    y <- 1 + 0.3 * x1 + 3 * x3 + e
    for (method in c("forward", "two.ways.forward")) {
        fit <- select_with(selectors[[method]], y, d)
        expect_identical(colnames(fit$data), c("x3", "x1"), info = method)
    }
    for (method in c("backward", "two.ways.backward")) {
        fit <- select_with(selectors[[method]], y, d)
        expect_identical(colnames(fit$data), c("x1", "x3"), info = method)
    }
})

test_that("alfa is compared as maSigPro does", {
    y <- 1 + 0.3 * x1 + 3 * x3 + e
    full <- summary(glm(y ~ ., data = d, epsilon = 1e-5))$coefficients
    single <- summary(glm(y ~ x3, data = d, epsilon = 1e-5))$coefficients
    # Backward keeps a term whose p-value equals alfa
    alfa <- full["x2", 4]
    expect_identical(
        names(coef(sc.stepback(y, d, alfa = alfa))),
        c("(Intercept)", "x1", "x2", "x3")
    )
    # Forward needs a p-value below alfa, two-ways forward at most alfa
    alfa <- single["x3", 4]
    expect_identical(names(coef(sc.stepfor(y, d, alfa = alfa))), "(Intercept)")
    expect_identical(
        names(coef(sc.two.ways.stepfor(y, d, alfa = alfa)))[1:2],
        c("(Intercept)", "x3")
    )
})

test_that("non-syntactic column names need no cleaning", {
    y <- 1 + 0.3 * x1 + 3 * x3 + e
    d_bt <- d
    colnames(d_bt) <- c("path A|time", "path B|x", "path C x")
    for (method in names(selectors)) {
        plain <- select_with(selectors[[method]], y, d)
        odd <- select_with(selectors[[method]], y, d_bt)
        expect_identical(unname(coef(odd)), unname(coef(plain)), info = method)
        expect_identical(
            colnames(odd$data),
            colnames(d_bt)[match(colnames(plain$data), colnames(d))],
            info = method
        )
        expect_identical(deviance(odd), deviance(plain), info = method)
    }
})

test_that("offset, weights, epsilon and max_it reach the fit", {
    y <- 1 + 2 * x1 + e
    off <- log(c(1, 2, 1, 2, 1, 2, 1, 2))
    w <- c(1, 2, 1, 2, 2, 1, 2, 1)
    for (method in names(selectors)) {
        fit <- selectors[[method]](
            y = y, d = d, alfa = 0.05, family = gaussian(),
            epsilon = 1e-6, useOffset = off, useWeight = w, max_it = 7
        )
        expect_identical(unname(fit$offset), off, info = method)
        expect_identical(unname(fit$prior.weights), w, info = method)
        expect_identical(fit$control$epsilon, 1e-6, info = method)
        expect_identical(fit$control$maxit, 7, info = method)
    }
})

test_that("a column that is not estimable is kept by backward only", {
    y <- 1 + 2 * x1 + e
    d_al <- cbind(d, copy = x1)
    expect_identical(
        names(coef(sc.stepback(y, d_al))),
        c("(Intercept)", "x1", "copy")
    )
    expect_true(is.na(coef(sc.stepback(y, d_al))[["copy"]]))
    for (method in c("forward", "two.ways.backward", "two.ways.forward")) {
        fit <- select_with(selectors[[method]], y, d_al)
        expect_identical(names(coef(fit)), c("(Intercept)", "x1"),
            info = method
        )
    }
})

test_that("backward methods keep the row names of the design", {
    y <- 1 + 0.3 * x1 + 3 * x3 + e
    expect_identical(rownames(sc.stepback(y, d)$model), rownames(d))
    expect_identical(rownames(sc.stepfor(y, d)$model), as.character(1:8))
    # With one column left, the rows are numbered (as in maSigPro)
    y <- 1 + 2 * x1 + e
    expect_identical(rownames(sc.stepback(y, d)$model), as.character(1:8))
})

test_that("extract_fitting() sets non-estimable terms to NA", {
    y <- 1 + 2 * x1 + e
    dis <- as.matrix(cbind(d, copy = x1))
    lmf <- glm(y ~ ., data = as.data.frame(dis))
    fit <- extract_fitting(
        reg = sc.stepback(y, as.data.frame(dis)), lmf = lmf,
        model.glm.0 = glm(y ~ 1), dis = dis, family = gaussian(),
        name = "g", vars.in = colnames(dis), alfa = 0.05,
        influ.info = matrix(NA, nrow = 8, ncol = 1)
    )
    expect_identical(is.na(fit$coeff), c(FALSE, FALSE, FALSE, FALSE, TRUE))
    expect_equal(fit$coeff[2], 2)
})

# Negative binomial counts on which the two-ways methods cycle or empty the
# selection (found by a random search, data rounded to two decimals)
nb <- MASS::negative.binomial(10)

test_that("two-ways forward selection stops when it cycles", {
    # x3 enters, then x4; x2 enters and makes x4 and x2 non-significant, so
    # both leave and x4 enters again. maSigPro's two.ways.stepfor() never
    # returns; the new code stops when a state comes back.
    y <- c(1, 4, 12, 10, 11, 9, 43, 6, 54, 2)
    d_cyc <- data.frame(
        x1 = c(
            0.52, -0.61, -0.44, -0.76, 0.53, -0.46, -0.35, -0.11, 0.05, -0.91
        ),
        x2 = c(
            0.96, -0.51, -0.36, -0.72, 0.84, -0.5, -0.94, 0.37, -0.59, -0.54
        ),
        x3 = c(
            -0.15, -0.14, -1.58, -0.49, 1.25, -0.2, 1.61, -2.64, 1.98, -0.83
        ),
        x4 = c(
            -2.53, -1.03, 0.11, -1.15, 0.36, -0.62, -0.4, -0.32, -0.63, -0.73
        )
    )
    fit <- sc.two.ways.stepfor(y, d_cyc, alfa = 0.1, family = nb)
    expect_identical(names(coef(fit)), c("(Intercept)", "x3", "x4"))
})

test_that("two-ways forward selection can remove every term", {
    # maSigPro's two.ways.stepfor() stops with an error on these data
    y <- c(3, 9, 41, 3, 1, 1, 5, 16, 5, 37, 2, 39)
    d_out <- data.frame(
        x1 = c(
            0.06, -0.37, -0.55, 1.39, -0.65, 0.04, 0.58, -0.29, 1.37, 0.61,
            -0.42, 0.39
        ),
        x2 = c(
            -1.68, 0.04, -0.09, -0.37, -1.35, -0.17, 1.56, 0.18, 0.66, 1.69,
            -1.23, 0.67
        ),
        x3 = c(
            -0.51, 0.42, 1.29, -0.78, 1.06, -0.01, -0.11, -0.04, 0.92, 0.78,
            0.63, 0.97
        )
    )
    fit <- sc.two.ways.stepfor(y, d_out, alfa = 0.05, family = nb)
    expect_identical(names(coef(fit)), "(Intercept)")
    expect_identical(
        names(coef(sc.stepfor(y, d_out, alfa = 0.05, family = nb))),
        c("(Intercept)", "x3", "x2")
    )
})

test_that("two-ways backward selection stops when it cycles", {
    # x3 (p just above 0.05) leaves, then x1; both come back and the cycle
    # starts again. maSigPro's two.ways.stepback() stops with an error here.
    y <- c(6, 1, 6, 3, 4, 11, 9, 50, 3, 11, 11, 7, 38, 13, 4, 25, 41, 18, 6, 2)
    d_cyc <- data.frame(
        x1 = c(
            -1.17, -1.16, 1.07, -0.26, -1.38, -0.35, -0.77, 0.32, 0.24, -0.47,
            0.17, 0.69, 0.96, 0.75, 2.74, 0.31, 0.2, -0.33, 1.05, 0.99
        ),
        x2 = c(
            -1.22, -0.87, 0.39, -0.73, -0.61, -0.1, -0.99, 0.37, 0.09, -1.16,
            0.57, 0.81, 0.42, 0.34, 2.87, 0.58, -0.36, -0.81, 1.05, 0.53
        ),
        x3 = c(
            -0.49, -1.23, -0.73, 0.88, 0.68, 0.49, -0.63, 0.96, 2.06, 0.09,
            -0.13, 0.52, 1.2, -0.49, 0.1, -1.53, -0.52, -0.37, -0.8, 1.12
        ),
        x4 = c(
            0.24, -1.04, -0.34, -0.84, -1.34, 0.37, -0.83, 0.77, -0.87, -0.44,
            0.79, 0.09, 0.17, -0.15, -0.39, 1.51, 0.74, 0.2, 0.15, -1.35
        )
    )
    fit <- sc.two.ways.stepback(y, d_cyc, alfa = 0.05, family = nb)
    expect_identical(names(coef(fit)), c("(Intercept)", "x2", "x4", "x3", "x1"))
})
