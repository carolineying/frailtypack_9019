suppressPackageStartupMessages({
  library(frailtypack)
  library(survival)
})

matrix_test <- getFromNamespace(".adaptiveGKMatrixTest", "frailtypack")
mode_test <- getFromNamespace(".adaptiveModeTest", "frailtypack")
survival_test <- getFromNamespace(".adaptiveSurvivalTest", "frailtypack")
ascertainment_test <- getFromNamespace(".ascertainmentGaussianTest",
                                       "frailtypack")

output_dir <- Sys.getenv("FRAILTYPACK_AUDIT_OUTPUT", unset = "audit-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add_result <- function(test, expected, observed, tolerance, pass,
                       abs_error = NA_real_, rel_error = NA_real_,
                       log_error = NA_real_, evaluations = NA_integer_,
                       details = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    test = test, expected = as.character(expected),
    observed = as.character(observed), tolerance = as.character(tolerance),
    absolute_error = abs_error, relative_error = rel_error,
    log_error = log_error, evaluations = evaluations,
    pass = isTRUE(pass), details = details, stringsAsFactors = FALSE
  )
}
check_scalar <- function(test, observed, expected, tolerance,
                         evaluations = NA_integer_, details = "") {
  ae <- abs(observed - expected)
  re <- ae / max(abs(expected), .Machine$double.xmin)
  le <- if (observed > 0 && expected > 0) abs(log(observed) - log(expected)) else Inf
  add_result(test, format(expected, digits = 17), format(observed, digits = 17),
             tolerance, ae <= tolerance, ae, re, le, evaluations, details)
}

spd <- function(d, scale = 1, rho = 0.25) {
  scale * ((1 - rho) * diag(d) + rho * outer(seq_len(d), seq_len(d),
                                              function(i, j) (-1)^(i + j) /
                                                sqrt(i * j)))
}

## Phases 4--5: normalization, including arbitrary mismatched proposals.
for (d in c(1L, 2L, 3L, 5L)) {
  priors <- list(identity = diag(d), correlated = spd(d, rho = 0.35))
  for (scale in c(0.5, 1, 2)) {
    for (prior_name in names(priors)) {
      prior <- scale * priors[[prior_name]]
      prior_chol <- t(chol(prior))
      proposals <- list(
        matched = list(mean = rep(0, d), cov = prior),
        shifted_scaled = list(mean = seq_len(d) / (20 * d),
                              cov = 1.2 * prior)
      )
      for (proposal_name in names(proposals)) {
        q <- proposals[[proposal_name]]
        ## The native correction is written in standardized prior coordinates.
        ## Convert the deliberately chosen b-space proposal to z-space.
        q_mean_z <- solve(prior_chol, q$mean)
        q_cov_z <- solve(prior_chol, q$cov)
        q_cov_z <- t(solve(prior_chol, t(q_cov_z)))
        a <- matrix_test(rep(0, d), diag(d), q_mean_z, q_cov_z,
                         test.case = "constant",
                         method = "adaptive-genz-keister",
                         maxpts = 7500L, epsrel = 1e-11)
        check_scalar(sprintf("constant adaptive d=%d prior=%s scale=%g q=%s",
                             d, prior_name, scale, proposal_name),
                     a$result, 1, if (d <= 3L) 2e-10 else 2e-7,
                     a$neval, sprintf("ifail=%d status=%d", a$ifail, a$status))
      }
      o <- matrix_test(rep(0, d), prior, test.case = "constant",
                       method = "genz-keister")
      check_scalar(sprintf("constant original d=%d prior=%s scale=%g",
                           d, prior_name, scale), o$result, 1, 2e-13, o$neval)
    }
  }
}

## A missing sqrt(2) changes this expectation by an order-one amount.
gk_scalar <- getFromNamespace(".adaptiveGKGaussianTest", "frailtypack")
gaussian_convention <- gk_scalar(1L, 2L, "genz-keister",
                                 shift = 0.8, ridge = 0.7)
gaussian_exact <- exp(0.5 * 0.8^2 / 1.7) / sqrt(1.7)
check_scalar("normalized-N(0,1) convention / sqrt(2) detector",
             gaussian_convention$result, gaussian_exact, 3e-11,
             gaussian_convention$neval)

## Phase 6: shifted non-diagonal Gaussian. With q equal to the target, the
## correction is constant and a one-node rule must return c.
for (d in 1:4) {
  mu <- (-1)^seq_len(d) * seq_len(d) / 4
  V <- spd(d, scale = 0.7, rho = 0.3)
  c0 <- 1.25 + d / 7
  a <- matrix_test(mu, V, mu, V, test.case = "shifted-gaussian",
                   constant = c0, minpts = 1L, maxpts = 1L,
                   epsabs = 1e-15, epsrel = 1e-15)
  check_scalar(sprintf("shifted Gaussian exact proposal d=%d one point", d),
               a$result, c0, 2e-13, a$neval,
               sprintf("ifail=%d status=%d", a$ifail, a$status))
}

## Helpers for independent direct and tensor-product references.
normal_gh <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- sqrt(seq_len(n - 1L))
  J <- J + t(J)
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1L, ]^2)
}
tensor_survival <- function(Sigma, delta, hazard, order = 20L) {
  d <- length(delta)
  gh <- normal_gh(order)
  grid <- expand.grid(rep(list(seq_along(gh$nodes)), d))
  z <- matrix(gh$nodes[as.matrix(grid)], ncol = d)
  weights <- matrix(gh$weights[as.matrix(grid)], ncol = d)
  b <- z %*% t(t(chol(Sigma)))
  value <- exp(drop(b %*% delta) - drop(exp(b) %*% hazard))
  sum(apply(weights, 1L, prod) * value)
}

## Phase 7: one-dimensional survival-integrand shapes.
one_d_cases <- list(
  approximately_gaussian = list(delta = 1, hazard = 1, variance = 0.5),
  shifted = list(delta = 3, hazard = 0.8, variance = 1),
  moderately_skewed = list(delta = 0, hazard = 0.35, variance = 1.5),
  broad_large_prior = list(delta = 1, hazard = 0.08, variance = 2)
)
one_d_table <- list()
for (nm in names(one_d_cases)) {
  p <- one_d_cases[[nm]]
  ref <- integrate(function(b)
    dnorm(b, sd = sqrt(p$variance)) *
      exp(p$delta * b - p$hazard * exp(b)),
    -Inf, Inf, rel.tol = 2e-13, subdivisions = 2000L)$value
  original <- survival_test(matrix(p$variance), p$delta, p$hazard,
                            method = "genz-keister", maxpts = 7500L,
                            epsrel = 1e-11)
  adaptive <- survival_test(matrix(p$variance), p$delta, p$hazard,
                            maxpts = 7500L, epsrel = 1e-11)
  one <- survival_test(matrix(p$variance), p$delta, p$hazard,
                       minpts = 1L, maxpts = 1L,
                       epsabs = 1e-15, epsrel = 1e-15)
  for (method_name in c("original", "adaptive", "laplace", "one-point")) {
    value <- switch(method_name, original = original$result,
                    adaptive = adaptive$result, laplace = adaptive$laplace,
                    `one-point` = one$result)
    evals <- switch(method_name, original = original$neval,
                    adaptive = adaptive$neval, laplace = 0L,
                    `one-point` = one$neval)
    one_d_table[[length(one_d_table) + 1L]] <- data.frame(
      shape = nm, method = method_name, reference = ref, estimate = value,
      absolute_error = abs(value - ref),
      relative_error = abs(value / ref - 1),
      log_error = abs(log(value) - log(ref)), evaluations = evals
    )
  }
  check_scalar(paste("one-point equals Laplace", nm),
               one$result, one$laplace, 5e-14, one$neval)
  check_scalar(paste("1D adaptive vs integrate", nm),
               adaptive$result, ref, 1e-9, adaptive$neval)
}
write.csv(do.call(rbind, one_d_table),
          file.path(output_dir, "one-dimensional-benchmark.csv"),
          row.names = FALSE)

## Phase 8: low-dimensional correlated and diagonal tensor references.
low_d_table <- list()
for (d in 2:4) {
  delta <- rep(c(1, 0), length.out = d)
  hazard <- seq(0.25, 0.8, length.out = d)
  for (cov_name in c("diagonal", "correlated")) {
    Sigma <- if (cov_name == "diagonal") diag(seq(0.6, 1.4, length.out = d)) else
      spd(d, scale = 1.2, rho = 0.3)
    ref <- tensor_survival(Sigma, delta, hazard, order = if (d == 4) 18L else 25L)
    original <- survival_test(Sigma, delta, hazard,
                              method = "genz-keister", maxpts = 7500L,
                              epsrel = 1e-10)
    adaptive <- survival_test(Sigma, delta, hazard, maxpts = 7500L,
                              epsrel = 1e-10)
    for (method_name in c("original", "adaptive", "laplace")) {
      value <- switch(method_name, original = original$result,
                      adaptive = adaptive$result, laplace = adaptive$laplace)
      evals <- switch(method_name, original = original$neval,
                      adaptive = adaptive$neval, laplace = 0L)
      low_d_table[[length(low_d_table) + 1L]] <- data.frame(
        dimension = d, covariance = cov_name, method = method_name,
        reference = ref, estimate = value, absolute_error = abs(value - ref),
        relative_error = abs(value / ref - 1),
        log_error = abs(log(value) - log(ref)), evaluations = evals,
        mode_gradient_norm = adaptive$gradnorm,
        hessian_min_eigenvalue = min(eigen(adaptive$hessian,
                                            symmetric = TRUE,
                                            only.values = TRUE)$values),
        hessian_condition = adaptive$condition,
        hessian_regularized = adaptive$regularized
      )
    }
    check_scalar(sprintf("tensor adaptive d=%d %s", d, cov_name),
                 adaptive$result, ref,
                 c(`2` = 5e-7, `3` = 3e-5, `4` = 6e-4)[as.character(d)],
                 adaptive$neval)
  }
}
write.csv(do.call(rbind, low_d_table),
          file.path(output_dir, "low-dimensional-benchmark.csv"),
          row.names = FALSE)

## Phase 9: proposal invariance on the same survival target.
d <- 3L
Sigma <- spd(d, scale = 1.3, rho = 0.3)
delta <- c(1, 0, 1)
hazard <- c(0.25, 0.7, 0.45)
proposal_specs <- list(
  q1 = list(scale = 1, shift = rep(0, d)),
  q2 = list(scale = sqrt(1.5), shift = rep(0, d)),
  q3 = list(scale = sqrt(0.75), shift = rep(0, d)),
  q4 = list(scale = 1, shift = c(0.12, -0.08, 0.05))
)
proposal_values <- vapply(proposal_specs, function(q)
  survival_test(Sigma, delta, hazard, proposal.scale = q$scale,
                proposal.shift = q$shift, maxpts = 7500L,
                epsrel = 1e-10)$result, numeric(1))
proposal_spread <- max(abs(proposal_values / proposal_values[[1L]] - 1))
add_result("proposal invariance survival d=3", "relative spread <= 1e-4",
           paste(format(proposal_values, digits = 12), collapse = ", "),
           "1e-4", proposal_spread <= 1e-4, rel_error = proposal_spread,
           details = "Loose survival tolerance reflects the 7,500-evaluation cap.")

## Phases 10--11: production mode, gradient, Hessian, Cholesky orientation.
d <- 4L
Sigma <- spd(d, scale = 1.1, rho = 0.32)
U <- t(chol(Sigma))
delta <- c(1, 0, 1, 0)
hazard <- c(0.2, 0.5, 0.8, 0.35)
m <- mode_test(U, delta, hazard)
hfun <- function(z) {
  b <- drop(U %*% z)
  sum(delta * b - hazard * exp(b)) - sum(z^2) / 2
}
fd_gradient <- function(f, x, eps = 2e-6)
  vapply(seq_along(x), function(i) {
    e <- rep(0, length(x)); e[i] <- eps
    (f(x + e) - f(x - e)) / (2 * eps)
  }, numeric(1))
fd_hessian <- function(f, x, eps = 2e-4) {
  d <- length(x); H <- matrix(0, d, d); f0 <- f(x)
  for (i in seq_len(d)) {
    ei <- rep(0, d); ei[i] <- eps
    H[i, i] <- (f(x + ei) - 2 * f0 + f(x - ei)) / eps^2
    if (i < d) for (j in (i + 1L):d) {
      ej <- rep(0, d); ej[j] <- eps
      H[i, j] <- H[j, i] <-
        (f(x + ei + ej) - f(x + ei - ej) -
           f(x - ei + ej) + f(x - ei - ej)) / (4 * eps^2)
    }
  }
  H
}
g_fd <- fd_gradient(hfun, m$mode.z)
H_fd_negative <- -fd_hessian(hfun, m$mode.z)
gradient_error <- max(abs(g_fd - m$gradient))
hessian_error <- max(abs(H_fd_negative - m$hessian))
symmetry_error <- max(abs(m$hessian - t(m$hessian)))
chol_error <- max(abs(m$hessian - m$chol.h %*% t(m$chol.h)))
eigenvalues <- eigen(m$hessian, symmetric = TRUE, only.values = TRUE)$values
set.seed(104729)
directions <- replicate(20L, {
  x <- rnorm(d); x / sqrt(sum(x^2))
})
local_drops <- apply(directions, 2L, function(v)
  hfun(m$mode.z) - hfun(m$mode.z + 1e-4 * v))
add_result("mode gradient finite difference", "max error <= 2e-7",
           format(gradient_error, digits = 8), "2e-7",
           gradient_error <= 2e-7, abs_error = gradient_error,
           details = sprintf("native gradnorm=%g iterations=%d status=%d",
                             m$gradnorm, m$iterations, m$status))
add_result("negative Hessian finite difference", "max error <= 2e-6",
           format(hessian_error, digits = 8), "2e-6",
           hessian_error <= 2e-6, abs_error = hessian_error,
           details = paste("eigenvalues", paste(signif(eigenvalues, 8),
                                                collapse = ",")))
add_result("Hessian symmetry", "max asymmetry <= 5e-15",
           format(symmetry_error, digits = 8), "5e-15",
           symmetry_error <= 5e-15, abs_error = symmetry_error)
add_result("Cholesky orientation H=C C'", "max error <= 5e-14",
           format(chol_error, digits = 8), "5e-14",
           chol_error <= 5e-14, abs_error = chol_error)
add_result("mode local perturbations decrease h", "all decreases > 0",
           format(min(local_drops), digits = 8), "> 0",
           all(local_drops > 0), details = paste("min/max drop",
                                                 signif(range(local_drops), 6),
                                                 collapse = "/"))

## Repeat the derivative comparison for several small family dimensions.
for (d_small in 1:3) {
  Sigma_small <- spd(d_small, scale = 0.8 + 0.1 * d_small, rho = 0.28)
  U_small <- t(chol(Sigma_small))
  delta_small <- rep(c(1, 0), length.out = d_small)
  hazard_small <- seq(0.25, 0.75, length.out = d_small)
  mode_small <- mode_test(U_small, delta_small, hazard_small)
  h_small <- function(z) {
    b <- drop(U_small %*% z)
    sum(delta_small * b - hazard_small * exp(b)) - sum(z^2) / 2
  }
  g_small_fd <- fd_gradient(h_small, mode_small$mode.z)
  H_small_fd <- -fd_hessian(h_small, mode_small$mode.z)
  g_small_error <- max(abs(g_small_fd - mode_small$gradient))
  H_small_error <- max(abs(H_small_fd - mode_small$hessian))
  add_result(sprintf("mode gradient finite difference d=%d", d_small),
             "max error <= 2e-7", format(g_small_error, digits = 8), "2e-7",
             g_small_error <= 2e-7, abs_error = g_small_error,
             details = sprintf("native gradnorm=%g status=%d",
                               mode_small$gradnorm, mode_small$status))
  add_result(sprintf("negative Hessian finite difference d=%d", d_small),
             "max error <= 2e-6", format(H_small_error, digits = 8), "2e-6",
             H_small_error <= 2e-6, abs_error = H_small_error)
}

## Phase 11b: covariance scale and source-level Cholesky orientation.
## The exact production transformations are b = U z and
## b = U (m + C^(-T) x), where H_z = C C'.  This test checks both the
## algebraic covariance and a Monte Carlo realization, and explicitly fails
## if either transformation behaves like an erroneous sqrt(2) scaling.
A_b <- U %*% solve(t(m$chol.h))
risk_at_mode <- hazard * exp(drop(U %*% m$mode.z))
H_b <- diag(risk_at_mode, d, d) + solve(Sigma)
adaptive_cov_expected <- solve(H_b)
adaptive_cov_transform <- A_b %*% t(A_b)
adaptive_cov_algebra_error <- max(abs(adaptive_cov_transform -
                                        adaptive_cov_expected))
add_result("adaptive proposal covariance A A'=H_b^-1",
           "max error <= 5e-13",
           format(adaptive_cov_algebra_error, digits = 8), "5e-13",
           adaptive_cov_algebra_error <= 5e-13,
           abs_error = adaptive_cov_algebra_error)

set.seed(8675309)
z_mc <- matrix(rnorm(200000L * d), ncol = d)
original_mc_cov <- cov(z_mc %*% t(U))
adaptive_mc_cov <- cov(z_mc %*% t(A_b))
relative_cov_error <- function(observed, expected)
  sqrt(sum((observed - expected)^2) / sum(expected^2))
original_cov_error <- relative_cov_error(original_mc_cov, Sigma)
adaptive_cov_error <- relative_cov_error(adaptive_mc_cov,
                                         adaptive_cov_expected)
original_doubled_error <- relative_cov_error(original_mc_cov, 2 * Sigma)
adaptive_doubled_error <- relative_cov_error(adaptive_mc_cov,
                                             2 * adaptive_cov_expected)
add_result("original empirical covariance scale",
           "relative error <= 0.01 and closer to Sigma than 2 Sigma",
           format(original_cov_error, digits = 8), "0.01",
           original_cov_error <= 0.01 &&
             original_cov_error < original_doubled_error,
           rel_error = original_cov_error,
           details = sprintf("relative error versus doubled covariance=%g",
                             original_doubled_error))
add_result("adaptive empirical covariance scale",
           "relative error <= 0.01 and closer to H_b^-1 than 2 H_b^-1",
           format(adaptive_cov_error, digits = 8), "0.01",
           adaptive_cov_error <= 0.01 &&
             adaptive_cov_error < adaptive_doubled_error,
           rel_error = adaptive_cov_error,
           details = sprintf("relative error versus doubled covariance=%g",
                             adaptive_doubled_error))

## Phase 12: permutation invariance of the complete correlated family target.
base <- survival_test(Sigma, delta, hazard, maxpts = 7500L, epsrel = 1e-10)
set.seed(90210)
for (i in 1:5) {
  p <- sample.int(d)
  permuted <- survival_test(Sigma[p, p], delta[p], hazard[p],
                            maxpts = 7500L, epsrel = 1e-10)
  check_scalar(sprintf("family permutation invariance %d", i),
               permuted$result, base$result, 1e-5, permuted$neval,
               "Numerical, not algebraic, tolerance under the capped d=4 product rule.")
}

## Phases 15--16: signed rules and Laplace identity.
d <- 5L
signed_mean <- seq(-0.4, 0.5, length.out = d)
signed_cov <- spd(d, scale = 0.45, rho = 0.2)
signed <- matrix_test(signed_mean, signed_cov,
                      proposal.mean = signed_mean,
                      proposal.cov = signed_cov,
                      test.case = "shifted-gaussian", constant = 1.7,
                      maxpts = 7500L, epsrel = 1e-12)
cancellation <- if (is.finite(signed$lognegative))
  exp(signed$lognegative - signed$logpositive) else 0
check_scalar("signed-weight shifted Gaussian d=5",
             signed$result, 1.7, 2e-5, signed$neval,
             sprintf("logpos=%g logneg=%g cancellation=%g status=%d",
                     signed$logpositive, signed$lognegative,
                     cancellation, signed$status))
one <- survival_test(matrix(1.7), 2, 0.4, minpts = 1L, maxpts = 1L,
                     epsabs = 1e-15, epsrel = 1e-15)
check_scalar("one-point adaptive/Laplace exact identity",
             one$result, one$laplace, 5e-14, one$neval)

## Phase 14: the separate proband marginal integral uses the same N(0,1)
## convention after b=sqrt(variance)*z.
for (variance in c(0.5, 1, 2)) {
  for (cumulative_hazard in c(0.1, 1, 5)) {
    observed <- ascertainment_test(cumulative_hazard, variance)
    expected <- integrate(function(z)
      dnorm(z) * exp(-cumulative_hazard * exp(sqrt(variance) * z)),
      -Inf, Inf, rel.tol = 2e-13)$value
    check_scalar(sprintf("ascertainment Gaussian variance=%g H=%g",
                         variance, cumulative_hazard),
                 observed, expected, 2e-7,
                 details = "Fixed 20-node legacy Gauss-Hermite tolerance.")
    add_result(sprintf("ascertainment probability range variance=%g H=%g",
                       variance, cumulative_hazard),
               "0 < p < 1", format(observed, digits = 17), "(0,1)",
               is.finite(observed) && observed > 0 && observed < 1)
  }
}

## Phase 18: deterministic repeatability.
repeat_values <- replicate(5L, survival_test(
  spd(3L, scale = 1.2, rho = 0.3), c(1, 0, 1), c(0.3, 0.7, 0.5),
  maxpts = 7500L, epsrel = 1e-10)$result)
add_result("deterministic repeatability", "bit-identical",
           paste(format(repeat_values, digits = 17), collapse = ", "),
           "identical()", length(unique(repeat_values)) == 1L)

results <- do.call(rbind, rows)
write.csv(results, file.path(output_dir, "audit-test-table.csv"),
          row.names = FALSE)
print(results, row.names = FALSE)
if (any(!results$pass)) {
  stop(sum(!results$pass), " rigorous adaptive quadrature audit test(s) failed")
}
