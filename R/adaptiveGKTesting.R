.adaptiveGKGaussianTest <- function(dimension, test.case = 1L,
                                    method = c("genz-keister",
                                               "adaptive-genz-keister",
                                               "genz-keister-fallback",
                                               "adaptive-genz-keister-staged",
                                               "adaptive-genz-keister-tempered"),
                                    shift = rep(0, dimension), ridge = 0,
                                    minpts = 30L, maxpts = 7500L,
                                    epsabs = 1e-12, epsrel = 1e-10) {
  method <- match.arg(method)
  dimension <- as.integer(dimension)
  stopifnot(length(dimension) == 1L, dimension >= 1L,
            length(shift) == dimension, ridge > -1,
            maxpts <= 7500L)
  z <- .Fortran(C_adaptive_gk_test,
                dimension,
                as.integer(test.case),
                as.integer(match(method, c("genz-keister",
                                           "adaptive-genz-keister",
                                           "genz-keister-fallback",
                                           "adaptive-genz-keister-staged",
                                           "adaptive-genz-keister-tempered")) - 1L),
                as.integer(minpts), as.integer(maxpts),
                as.double(epsabs), as.double(epsrel),
                as.double(shift), as.double(ridge),
                result = double(1), abserr = double(1),
                neval = integer(1), ifail = integer(1),
                analytic = double(1), status = integer(1))
  unclass(z[c("result", "analytic", "abserr", "neval", "ifail", "status")])
}

.correlatedWeibullCleanupTest <- function() {
  unname(.Fortran(C_correlated_weib_cleanup_test, status = integer(1L))$status)
}

.adaptiveGKMatrixTest <- function(target.mean, target.cov,
                                  proposal.mean = target.mean,
                                  proposal.cov = target.cov,
                                  proposal.scale = 1,
                                  test.case = c("shifted-gaussian", "constant"),
                                  method = c("adaptive-genz-keister",
                                             "genz-keister"),
                                  constant = 1, minpts = 30L,
                                  maxpts = 7500L, epsabs = 1e-12,
                                  epsrel = 1e-10) {
  test.case <- match.arg(test.case)
  method <- match.arg(method)
  d <- length(target.mean)
  stopifnot(d >= 1L, identical(dim(target.cov), c(d, d)),
            length(proposal.mean) == d,
            identical(dim(proposal.cov), c(d, d)),
            proposal.scale > 0, constant > 0, maxpts <= 7500L)
  target.precision.chol <- t(chol(solve(target.cov)))
  proposal.precision.chol <- t(chol(solve(proposal.cov)))
  z <- .Fortran(
    C_adaptive_gk_matrix_test,
    as.integer(d),
    as.integer(match(test.case, c("constant", "shifted-gaussian"))),
    as.integer(match(method, c("genz-keister",
                               "adaptive-genz-keister")) - 1L),
    as.integer(minpts), as.integer(maxpts),
    as.double(epsabs), as.double(epsrel),
    as.double(target.mean), as.double(target.precision.chol),
    as.double(proposal.mean), as.double(proposal.precision.chol),
    as.double(proposal.scale), as.double(constant),
    result = double(1), abserr = double(1), neval = integer(1),
    ifail = integer(1), status = integer(1),
    logpositive = double(1), lognegative = double(1)
  )
  unclass(z[c("result", "abserr", "neval", "ifail", "status",
              "logpositive", "lognegative")])
}

.adaptiveModeTest <- function(prior.chol, delta, hazard, mode.tol = 1e-10,
                              maxit = 100L, hessian.eps = 1e-10) {
  d <- length(delta)
  stopifnot(identical(dim(prior.chol), c(d, d)),
            length(hazard) == d, all(hazard >= 0))
  z <- .Fortran(
    C_adaptive_mode_test, as.integer(d), as.double(prior.chol),
    as.double(delta), as.double(hazard), as.double(mode.tol),
    as.integer(maxit), as.double(hessian.eps),
    mode.z = double(d), log.mode = double(1), gradient = double(d),
    hessian = double(d * d), chol.h = double(d * d),
    status = integer(1), iterations = integer(1), gradnorm = double(1),
    mineig = double(1), minstab = double(1), condition = double(1),
    regularized = integer(1)
  )
  z$hessian <- matrix(z$hessian, d, d)
  z$chol.h <- matrix(z$chol.h, d, d)
  unclass(z[c("mode.z", "log.mode", "gradient", "hessian", "chol.h",
              "status", "iterations", "gradnorm", "mineig", "minstab",
              "condition", "regularized")])
}

.adaptiveSurvivalTest <- function(prior.cov, delta, hazard,
                                  method = c("adaptive-genz-keister",
                                             "genz-keister"),
                                  proposal.scale = 1,
                                  proposal.shift = rep(0, length(delta)),
                                  minpts = 30L, maxpts = 7500L,
                                  epsabs = 1e-12, epsrel = 1e-10) {
  method <- match.arg(method)
  d <- length(delta)
  stopifnot(identical(dim(prior.cov), c(d, d)),
            length(hazard) == d, all(hazard >= 0),
            length(proposal.shift) == d, proposal.scale > 0,
            maxpts <= 7500L)
  prior.chol <- t(chol(prior.cov))
  z <- .Fortran(
    C_adaptive_survival_test, as.integer(d),
    as.integer(match(method, c("genz-keister",
                               "adaptive-genz-keister")) - 1L),
    as.integer(minpts), as.integer(maxpts), as.double(epsabs),
    as.double(epsrel), as.double(prior.chol), as.double(delta),
    as.double(hazard), as.double(proposal.scale),
    as.double(proposal.shift), result = double(1), abserr = double(1),
    neval = integer(1), ifail = integer(1), status = integer(1),
    laplace = double(1), mode.z = double(d), gradient = double(d),
    hessian = double(d * d), chol.h = double(d * d),
    iterations = integer(1), gradnorm = double(1), mineig = double(1),
    minstab = double(1), condition = double(1), regularized = integer(1)
  )
  z$hessian <- matrix(z$hessian, d, d)
  z$chol.h <- matrix(z$chol.h, d, d)
  unclass(z[c("result", "abserr", "neval", "ifail", "status", "laplace",
              "mode.z", "gradient", "hessian", "chol.h", "iterations",
              "gradnorm", "mineig", "minstab", "condition", "regularized")])
}

.ascertainmentGaussianTest <- function(cumulative.hazard, variance) {
  stopifnot(length(cumulative.hazard) == 1L, cumulative.hazard >= 0,
            length(variance) == 1L, variance >= 0)
  unname(.Fortran(C_ascertainment_gaussian_test,
                  as.double(cumulative.hazard), as.double(variance),
                  result = double(1))$result)
}
