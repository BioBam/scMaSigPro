#' Estimated theta for the negative binomial family
#'
#' Internal helpers behind `sc.p.vector(theta_estimate = "gene")` and
#' `sc.p.vector(theta_estimate = "global")` (decision D16). By default the
#' negative binomial family has a fixed theta (`negative.binomial(10)`), as
#' in maSigPro. On request, theta is estimated by maximum likelihood on each
#' gene's full model: the fit and theta are updated in turn (as
#' [MASS::glm.nb()] does) until theta changes by less than 1e-6 on the log
#' scale, at most `max_it` times. Each theta step maximises the negative
#' binomial log-likelihood with the fitted means held fixed, within
#' `.theta_bounds`; a gene whose counts are no more variable than Poisson
#' counts gets the upper bound (a Poisson-like fit) instead of the error or
#' the endless growth of [MASS::theta.ml()]. With `"gene"` each gene is then
#' fitted, full and intercept-only model alike, with its own theta; with
#' `"global"` every gene is fitted with the median of the estimates.
#'
#' @importFrom stats optimize median glm
#' @importFrom MASS negative.binomial
#' @name theta-estimate
#' @noRd
NULL

# Search interval for theta: 1e-2 is extreme overdispersion, 1e6 is
# Poisson-like for counts up to about 1e5
.theta_bounds <- c(1e-2, 1e6)

# TRUE for the families of MASS::negative.binomial()
.is_nb_family <- function(family) {
    is.character(family$family) && length(family$family) == 1L &&
        startsWith(family$family, "Negative Binomial(")
}

# The family `family` (a negative binomial family) with theta set to
# `theta`; the link functions and the link label are those of `family`
.nb_with_theta <- function(family, theta) {
    fam <- negative.binomial(theta)
    link_parts <- c("link", "linkfun", "linkinv", "mu.eta", "valideta")
    fam[link_parts] <- family[link_parts]
    fam
}

# The family of one gene: its own theta when sc.p.vector() estimated one
# theta per gene, otherwise `family` itself
.gene_family <- function(family, gene) {
    if (identical(family$theta_estimate, "gene")) {
        return(.nb_with_theta(family, family$theta_gene[[gene]]))
    }
    family
}

# Maximum likelihood theta of gene i (row i of dat) on the full model. The
# start is the theta of `family`. Returns the estimate, whether it is at a
# bound of .theta_bounds, whether the alternation converged, and the number
# of alternations. sc.p.vector() runs it through bplapply(), so it takes all
# its data as arguments: SnowParam workers are new R sessions.
.theta_gene <- function(i, dat, dis, family, epsilon, offset_data, max_it) {
    y <- as.numeric(dat[i, ])
    # The bins with a value: glm() leaves the others out of the fit, so the
    # likelihood pairs these with the fitted means
    y_obs <- y[!is.na(y)]
    loglik <- function(log_theta, y, mu) {
        th <- exp(log_theta)
        ll <- lgamma(y + th) - lgamma(th) - lgamma(y + 1) + th * log(th) +
            y * log(mu) - (y + th) * log(th + mu)
        sum(ll)
    }
    theta <- environment(family$variance)$.Theta
    converged <- FALSE
    n_alt <- 0L
    for (k in seq_len(max_it)) {
        n_alt <- k
        fit <- glm(y ~ .,
            data = dis, family = .nb_with_theta(family, theta),
            epsilon = epsilon, offset = offset_data, maxit = max_it
        )
        opt <- optimize(loglik, log(.theta_bounds),
            y = y_obs, mu = fit$fitted.values, maximum = TRUE, tol = 1e-10
        )
        new <- exp(opt$maximum)
        step <- abs(log(new / theta))
        theta <- new
        if (step < 1e-6) {
            converged <- TRUE
            break
        }
    }
    at_bound <- theta > .theta_bounds[2] / 1.001 ||
        theta < .theta_bounds[1] * 1.001
    c(theta = theta, at_bound = at_bound, converged = converged, n_alt = n_alt)
}
