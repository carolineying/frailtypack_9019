suppressPackageStartupMessages({
  library(frailtypack)
  library(survival)
})

gk_test <- getFromNamespace(".adaptiveGKGaussianTest", "frailtypack")
cleanup_test <- getFromNamespace(".correlatedWeibullCleanupTest", "frailtypack")
stopifnot(identical(cleanup_test(), 0L), identical(cleanup_test(), 0L))

## 1. Gaussian identity: the normalized Gaussian weight integrates one.
for (d in c(1L, 5L, 8L)) {
  for (method in c("genz-keister", "adaptive-genz-keister")) {
    z <- gk_test(d, test.case = 1L, method = method)
    stopifnot(abs(z$result - 1) < 2e-13, z$status == 0L)
  }
}

## 2. Shifted Gaussian product with a known normalizing constant.
for (d in c(1L, 2L, 4L, 5L, 8L)) {
  shift <- seq_len(d) / 5
  z <- gk_test(d, test.case = 2L, method = "adaptive-genz-keister",
               shift = shift, ridge = 0.7)
  stopifnot(abs(z$result / z$analytic - 1) < 2e-12,
            z$neval <= 7500L, z$status == 0L)
}

## 3. One-dimensional high-accuracy direct integration.
a <- 0.8
r <- 0.7
direct <- integrate(function(x) dnorm(x) * exp(a * x - r * x^2 / 2),
                    -Inf, Inf, rel.tol = 1e-13)$value
one.original <- gk_test(1L, 2L, "genz-keister", shift = a, ridge = r)
one.adaptive <- gk_test(1L, 2L, "adaptive-genz-keister",
                        shift = a, ridge = r)
stopifnot(abs(one.original$result - direct) < 2e-10,
          abs(one.adaptive$result - direct) < 2e-12)

## 4. Independent tensor-product Gauss--Hermite check in dimensions 2--4.
normal_gh <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- sqrt(seq_len(n - 1L))
  J <- J + t(J)
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1L, ]^2)
}
gh <- normal_gh(20L)
for (d in 2:4) {
  shift <- seq_len(d) / 7
  grid <- expand.grid(rep(list(seq_along(gh$nodes)), d))
  node_matrix <- matrix(gh$nodes[as.matrix(grid)], ncol = d)
  weight_matrix <- matrix(gh$weights[as.matrix(grid)], ncol = d)
  tensor <- sum(apply(weight_matrix, 1L, prod) *
                  exp(drop(node_matrix %*% shift) -
                        0.35 * rowSums(node_matrix^2)))
  z0 <- gk_test(d, 2L, "genz-keister", shift = shift, ridge = 0.7)
  z1 <- gk_test(d, 2L, "adaptive-genz-keister", shift = shift, ridge = 0.7)
  stopifnot(abs(tensor / z0$analytic - 1) < 2e-11,
            abs(z1$result / tensor - 1) < 2e-11)
}

## A deterministic correlated-Weibull fixture for interface, diagnostics,
## backward compatibility, and ascertainment smoke tests.
fixture <- data.frame(
  t0 = 0,
  time = c(1, 2, 3, 4, 1.2, 2.2, 3.2, 4.2, 1.4, 2.4, 3.4, 4.4),
  status = rep(c(1, 0, 1, 0), 3),
  x = rep(c(-0.5, 0.2, 0.7, -0.1), 3),
  fam = rep(1:3, each = 4)
)
K <- kronecker(diag(3), matrix(0.25, 4, 4))
diag(K) <- 1
base_args <- list(
  formula = Surv(t0, time, status) ~ x + cluster(fam),
  data = fixture, hazard = "Weibull", RandDist = "LogN",
  recurrentAG = TRUE, covMatrix1 = K, print.times = FALSE,
  maxit = 1L, init.B = 0,
  cubature.control = list(maxpts = 200L, epsabs = 1e-8, epsrel = 1e-6)
)

## 5. The omitted method and explicit legacy method are bit-for-bit equal.
legacy.default <- suppressWarnings(do.call(frailtyPenal, base_args))
legacy.explicit <- suppressWarnings(do.call(
  frailtyPenal, c(base_args, list(integration.method = "genz-keister"))))
stopifnot(identical(legacy.default$b, legacy.explicit$b),
          identical(legacy.default$logLik, legacy.explicit$logLik),
          identical(legacy.explicit$integration.failure, "fail"),
          legacy.explicit$integration.summary$floor.applications == 0L,
          legacy.explicit$integration.floor$applications == 0L,
          nrow(legacy.explicit$integration.floor$by.family) == 0L,
          isTRUE(legacy.explicit$integration.valid),
          identical(legacy.explicit$integration.status.label,
                    "budget-exhausted"))
bad.failure.action <- try(do.call(
  frailtyPenal, c(base_args, list(integration.failure = "not-an-action"))),
  silent = TRUE)
stopifnot(inherits(bad.failure.action, "try-error"))

## 6. Adaptive family diagnostics are finite and honor the hard budget.
adaptive <- suppressWarnings(do.call(
  frailtyPenal,
  c(base_args, list(integration.method = "adaptive-genz-keister",
                    adaptive.control = list(diagnostics = TRUE)))))
d <- adaptive$family.integration.diagnostics
stopifnot(all(c("covariance.factorization.status", "signed.sum.status",
                "signed.log.positive", "signed.log.negative",
                "signed.cancellation.ratio", "floor.applied", "floor.count",
                "floor.signed.nonpositive.count",
                "floor.nonpositive.value.count",
                "floor.prior.under.threshold.count",
                "floor.adaptive.under.threshold.count") %in% names(d)),
          nrow(d) == 3L, all(d$mode.status == 1L),
          all(is.finite(d$log.integral)), all(d$evaluations <= 200L),
          all(d$covariance.factorization.status == 0L),
          adaptive$integration.summary$maximum.dimension == 4L,
          adaptive$integration.summary$floor.applications == 0L,
          adaptive$integration.summary$families.ever.floored == 0L,
          all(d$floor.count == 0L), !any(d$floor.applied),
          adaptive$integration.summary$covariance.factorization.failures == 0L,
          isTRUE(adaptive$integration.valid),
          adaptive$integration.status %in% 0:1,
          !any(d$hard.failure))

## 6b. A forced mode failure retains its last finite Newton proposal.  This is
## still an exact quadrature coordinate change, and the reported evaluations
## never exceed maxpts.
nf <- 2L
m <- 12L
n <- nf * m
fallback.fixture <- data.frame(
  t0 = 0,
  time = rep(seq(1, 4, length.out = m), nf) + rep(c(0, 0.1), each = m),
  status = rep(c(1, 0), n / 2L),
  x = seq(-0.5, 0.7, length.out = n),
  fam = rep(seq_len(nf), each = m)
)
B <- matrix(0.25, m, m)
diag(B) <- 1
fallback.K <- kronecker(diag(nf), B)
forced.fallback <- suppressWarnings(frailtyPenal(
  Surv(t0, time, status) ~ x + cluster(fam),
  data = fallback.fixture, hazard = "Weibull", RandDist = "LogN",
  recurrentAG = TRUE, covMatrix1 = fallback.K, print.times = FALSE,
  maxit = 1L, init.B = 0,
  integration.method = "adaptive-genz-keister",
  cubature.control = list(maxpts = 200L, epsabs = 1e-100,
                          epsrel = 1e-100),
  adaptive.control = list(mode.tol = 1e-16, mode.maxit = 1L,
                          fallback = TRUE, diagnostics = TRUE,
                          warm.start = FALSE)
))
fd <- forced.fallback$family.integration.diagnostics
stopifnot(all(fd$mode.status < 0L), all(fd$fallback.used),
          all(fd$method == "adaptive-approximate-mode"),
          all(fd$evaluations <= 200L), all(fd$finite.integral),
          all(fd$budget.exhausted), !any(fd$hard.failure),
          isTRUE(forced.fallback$integration.valid),
          identical(forced.fallback$integration.status.label,
                    "fallback-used"))

## 6c. Nonfinite legacy/fallback values are rejected rather than floored to a
## small positive integral.  The test harness calls the production validator.
nonfinite <- suppressWarnings(gk_test(
  1L, 2L, "genz-keister", shift = 1000, ridge = 0,
  minpts = 30L, maxpts = 200L, epsabs = 1e-100, epsrel = 1e-100
))
stopifnot(nonfinite$status == -3L,
          !is.finite(nonfinite$result) || !is.finite(nonfinite$abserr))

## 6d. The same prior-centred rule is finite when invoked as a fallback,
## because its likelihood contribution and rule sum stay on the log scale.
## The analytic integral remains below double overflow, while individual raw
## node values overflow in the historical direct branch.
fallback.stable <- suppressWarnings(gk_test(
  1L, 2L, "genz-keister-fallback", shift = 169, ridge = 20,
  minpts = 30L, maxpts = 7500L, epsabs = 1e-100, epsrel = 1e-8
))
direct.unstable <- suppressWarnings(gk_test(
  1L, 2L, "genz-keister", shift = 169, ridge = 20,
  minpts = 30L, maxpts = 7500L, epsabs = 1e-100, epsrel = 1e-8
))
stopifnot(is.finite(fallback.stable$result), fallback.stable$result > 0,
          fallback.stable$status == 0L,
          direct.unstable$status == -3L)

## 6e. When both branches complete the same fully symmetric rule, changing
## only the arithmetic to signed log summation preserves the integral.
for (dimension in c(5L, 8L)) {
  shift <- seq_len(dimension) / 7
  direct.same.rule <- gk_test(
    dimension, 2L, "genz-keister", shift = shift, ridge = 0.7,
    minpts = 30L, maxpts = 7500L, epsabs = 1e-100, epsrel = 0.02
  )
  fallback.same.rule <- gk_test(
    dimension, 2L, "genz-keister-fallback", shift = shift, ridge = 0.7,
    minpts = 30L, maxpts = 7500L, epsabs = 1e-100, epsrel = 0.02
  )
  stopifnot(direct.same.rule$neval == fallback.same.rule$neval,
            abs(direct.same.rule$result - fallback.same.rule$result) /
              fallback.same.rule$analytic < 1e-12)
}

## 6f. Splitting a stable adaptive call into a staged call and restart must
## reproduce the same completed rule, while a broadened proposal remains an
## exact change of variables.
for (dimension in c(5L, 8L, 12L)) {
  shift <- seq_len(dimension) / 7
  single <- gk_test(
    dimension, 2L, "adaptive-genz-keister", shift = shift, ridge = 0.7,
    minpts = 30L, maxpts = 7500L, epsabs = 1e-100, epsrel = 1e-100
  )
  staged <- gk_test(
    dimension, 2L, "adaptive-genz-keister-staged", shift = shift,
    ridge = 0.7, minpts = 30L, maxpts = 7500L, epsabs = 1e-100,
    epsrel = 1e-100
  )
  stopifnot(single$neval == staged$neval,
            abs(single$result - staged$result) / single$analytic < 1e-12,
            staged$neval <= 7500L)
}
for (dimension in 1:2) {
  shift <- seq_len(dimension) / 7
  tempered <- gk_test(
    dimension, 2L, "adaptive-genz-keister-tempered", shift = shift,
    ridge = 0.7, minpts = 30L, maxpts = 7500L, epsabs = 1e-100,
    epsrel = 1e-8
  )
  stopifnot(abs(tempered$result / tempered$analytic - 1) < 1e-7,
            tempered$neval <= 7500L)
}

## 7. Ascertainment correction remains on its separate probability path and
## stays finite before the logarithm after the existing [eps,1-eps] safeguard.
asc_args <- base_args
## One proband in every family.
asc_args$proband <- as.integer(ave(fixture$fam, fixture$fam,
                                   FUN = seq_along) == 1L)
asc_args$currentage <- fixture$time + 0.5
asc <- suppressWarnings(do.call(
  frailtyPenal,
  c(asc_args, list(integration.method = "adaptive-genz-keister"))))
stopifnot(all(is.finite(asc$b)), !is.null(asc$integration.summary))

## corrRE=3: one covariance matrix whose variance is fixed.  This must use
## the same initialized proband-marginal ascertainment integral as corrRE=1.
fixed_asc_args <- asc_args
fixed_asc_args$sig.fixed <- 0.75
fixed_asc_1 <- suppressWarnings(do.call(
  frailtyPenal,
  c(fixed_asc_args, list(integration.method = "adaptive-genz-keister"))))
fixed_asc_2 <- suppressWarnings(do.call(
  frailtyPenal,
  c(fixed_asc_args, list(integration.method = "adaptive-genz-keister"))))
stopifnot(all(is.finite(fixed_asc_1$b)),
          is.finite(fixed_asc_1$logLik),
          identical(fixed_asc_1$b, fixed_asc_2$b),
          identical(fixed_asc_1$logLik, fixed_asc_2$logLik),
          !is.null(fixed_asc_1$integration.summary))

## A singular family covariance must be rejected before a partial Cholesky
## factor is passed to the conditional mode finder or quadrature callback.
singular_args <- base_args
singular_args$covMatrix1 <- kronecker(diag(3), matrix(1, 4, 4))
singular <- suppressWarnings(do.call(
  frailtyPenal,
  c(singular_args, list(
    integration.method = "adaptive-genz-keister",
    adaptive.control = list(diagnostics = TRUE)))))
stopifnot(!isTRUE(singular$integration.valid),
          singular$integration.summary$covariance.factorization.failures >= 1L,
          any(singular$family.integration.diagnostics$
                covariance.factorization.status != 0L),
          any(singular$family.integration.diagnostics$signed.sum.status == -9L),
          any(singular$family.integration.diagnostics$hard.failure))

## The two-covariance Weibull path sums its covariance components before the
## same adaptive family integral.
two_cov_args <- base_args
two_cov_args$covMatrix2 <- diag(nrow(fixture))
two_cov_args$sig.fixed <- 0.25
two_cov <- suppressWarnings(do.call(
  frailtyPenal,
  c(two_cov_args, list(integration.method = "adaptive-genz-keister",
                       adaptive.control = list(diagnostics = TRUE)))))
stopifnot(all(two_cov$family.integration.diagnostics$mode.status == 1L),
          all(is.finite(two_cov$family.integration.diagnostics$log.integral)))
