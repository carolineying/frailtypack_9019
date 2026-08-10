#!/usr/bin/env Rscript

## Small paired diagnostic driver. It intentionally refuses budgets above 7,500
## and does not launch a full simulation campaign.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 6L) {
  stop(paste(
    "Usage: Rscript compare_adaptive_genz_keister.R",
    "INPUT.rds OUTPUT.tsv [NFAM=10] [MAXIT=35] [RELATIONS=all] [PERMUTE=false]"
  ))
}
input <- args[[1L]]
output <- args[[2L]]
nfam <- if (length(args) >= 3L) as.integer(args[[3L]]) else 10L
maxit <- if (length(args) >= 4L) as.integer(args[[4L]]) else 35L
relations <- if (length(args) >= 5L) args[[5L]] else "all"
permute <- length(args) >= 6L && tolower(args[[6L]]) %in% c("true", "t", "1", "yes")
stopifnot(file.exists(input), nfam >= 1L, maxit >= 1L)

suppressPackageStartupMessages({
  library(frailtypack)
  library(survival)
})

z <- readRDS(input)
family <- z$family
## `simfam(..., variation = "kinship")` stores the conventional kinship
## matrix whose diagonal is 1/2.  The correlated-frailty simulation drivers
## fit `2 * kmat` so that the diagonal scale equals the reported frailty
## variance.  Newer corrected fixtures may store that model matrix directly
## as `K`.
K <- if (!is.null(z$K)) z$K else 2 * as.matrix(z$kmat)
if (relations != "all") {
  wanted <- as.integer(strsplit(relations, ",", fixed = TRUE)[[1L]])
  keep <- family$relation %in% wanted
  family <- family[keep, , drop = FALSE]
  K <- K[keep, keep, drop = FALSE]
}
ids <- unique(family$famID)[seq_len(min(nfam, length(unique(family$famID))))]
keep <- family$famID %in% ids
family <- family[keep, , drop = FALSE]
K <- K[keep, keep, drop = FALSE]
if (permute) {
  set.seed(73019L)
  p <- sample.int(nrow(family))
  family <- family[p, , drop = FALSE]
  K <- K[p, p, drop = FALSE]
}
family$t0 <- if ("t0" %in% names(family)) family$t0 else 0

gender_name <- if ("gender_cov" %in% names(family)) "gender_cov" else "gender"
gene_name <- if ("mgene_cov" %in% names(family)) "mgene_cov" else "mgene"
model_formula <- as.formula(sprintf(
  "Surv(t0, time, status) ~ %s + %s + cluster(famID)",
  gender_name, gene_name
))
true_sigma2 <- if (!is.null(z$settings$tau_k)) z$settings$tau_k else
  if (!is.null(z$true_sigma2)) z$true_sigma2 else 2

fit_one <- function(method, maxpts) {
  warnings <- character()
  elapsed <- system.time({
    fit <- withCallingHandlers(
      frailtyPenal(
        model_formula, data = family, hazard = "Weibull", RandDist = "LogN",
        print.times = FALSE, init.B = c(1, 0.5),
        init.Theta = sqrt(true_sigma2), covMatrix1 = K,
        recurrentAG = TRUE, maxit = maxit,
        proband = if ("proband" %in% names(family)) family$proband else NULL,
        currentage = if ("currentage" %in% names(family)) family$currentage else NULL,
        integration.method = method,
        adaptive.control = list(mode.tol = 1e-7, mode.maxit = 80L,
                                hessian.eps = 1e-10, fallback = TRUE,
                                diagnostics = TRUE, warm.start = TRUE),
        cubature.control = list(minpts = 30L, maxpts = maxpts,
                                epsabs = 1e-100, epsrel = 1e-6)
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })
  s <- fit$integration.summary
  outer_converged <- identical(as.integer(fit$istop), 1L)
  integration_valid <- isTRUE(fit$integration.valid)
  estimate_valid <- outer_converged && integration_valid
  sigma2_raw <- as.numeric(fit$sigma2)
  sigma2_hat <- if (estimate_valid) sigma2_raw else NA_real_
  max_evaluations <- if (is.null(s)) {
    unname(fit$cubature.diagnostics[["max_neval"]])
  } else {
    max(fit$family.integration.diagnostics$evaluations)
  }
  if (is.finite(max_evaluations) && max_evaluations > 7500L) {
    stop("A family used ", max_evaluations,
         " evaluations, exceeding the absolute 7,500-evaluation limit.")
  }
  data.frame(
    replicate = if (!is.null(z$replicate)) z$replicate else basename(input),
    permuted = permute,
    method = method, nominal_maxpts = maxpts,
    true_sigma2 = true_sigma2, sigma2_hat = sigma2_hat,
    sigma2_raw = sigma2_raw,
    bias = sigma2_hat - true_sigma2,
    squared_error = (sigma2_hat - true_sigma2)^2,
    converged = outer_converged, result_valid = estimate_valid,
    integration_failure = !integration_valid,
    integration_status = fit$integration.status.label,
    adaptive_completed = isTRUE(fit$adaptive.integration.completed),
    reached_cap = if (is.null(s)) NA_integer_ else s$reached.evaluation.cap,
    budget_exhausted = if (is.null(s)) NA_integer_ else s$budget.exhausted,
    fallback_families = if (is.null(s)) NA_integer_ else s$fallback.families,
    hard_failures = if (is.null(s)) NA_integer_ else s$hard.failures,
    nonfinite_integrals = if (is.null(s)) NA_integer_ else s$nonfinite.integrals,
    max_evaluations = max_evaluations,
    mode_failures = if (is.null(s)) NA_integer_ else s$mode.failures,
    hessian_regularizations = if (is.null(s)) NA_integer_ else s$hessian.regularizations,
    runtime_seconds = unname(elapsed[["elapsed"]]),
    log_likelihood = as.numeric(fit$logLik),
    fixed_effects = paste(
      paste(names(unlist(fit$coef)), format(unlist(fit$coef), digits = 12),
            sep = "="),
      collapse = ";"
    ),
    warning_count = length(unique(warnings)),
    n = nrow(family), families = length(ids),
    max_family_size = max(table(family$famID)),
    stringsAsFactors = FALSE
  )
}

## The historical direct method can start an indivisible rule that overruns a
## nominal 7,500 cap.  Keep its package-default 750-point comparison (whose
## completed rule remains below 7,500 here) and compare adaptive coordinates
## at both 750 and the strictly enforced 7,500 cap.
spec <- data.frame(
  method = c("genz-keister", "adaptive-genz-keister",
             "adaptive-genz-keister"),
  maxpts = c(750L, 750L, 7500L), stringsAsFactors = FALSE
)
rows <- lapply(seq_len(nrow(spec)), function(i) {
  fit_one(spec$method[[i]], spec$maxpts[[i]])
})

## External Laplace benchmark. This is not substituted into frailtypack's
## Weibull/ascertainment likelihood; it is a separate comparison only.
if (requireNamespace("coxme", quietly = TRUE)) {
  family$.row_id <- factor(seq_len(nrow(family)))
  dimnames(K) <- list(levels(family$.row_id), levels(family$.row_id))
  fixed_terms <- c(gender_name, gene_name)
  fixed_terms <- fixed_terms[vapply(family[fixed_terms], function(v)
    length(unique(v[!is.na(v)])) > 1L, logical(1))]
  cox_formula <- as.formula(sprintf(
    "Surv(t0, time, status) ~ %s(1 | .row_id)",
    if (length(fixed_terms)) paste0(paste(fixed_terms, collapse = " + "), " + ") else ""
  ))
  cox_error <- NULL
  elapsed <- system.time({
    cfit <- tryCatch(
      coxme::coxme(cox_formula, data = family,
                   varlist = coxme::coxmeMlist(list(K = K), rescale = FALSE)),
      error = function(e) { cox_error <<- conditionMessage(e); NULL }
    )
  })
  estimate <- if (is.null(cfit)) NA_real_ else
    as.numeric(unlist(cfit$vcoef, use.names = FALSE)[1L])
  cox_fixed <- if (is.null(cfit)) NA_character_ else
    paste(paste(names(coxme::fixef(cfit)),
                format(coxme::fixef(cfit), digits = 12), sep = "="),
          collapse = ";")
  rows[[length(rows) + 1L]] <- data.frame(
    replicate = if (!is.null(z$replicate)) z$replicate else basename(input),
    permuted = permute,
    method = "coxme-laplace", nominal_maxpts = NA_integer_,
    true_sigma2 = true_sigma2, sigma2_hat = estimate,
    bias = estimate - true_sigma2, squared_error = (estimate - true_sigma2)^2,
    sigma2_raw = estimate,
    converged = !is.null(cfit), result_valid = !is.null(cfit),
    integration_failure = NA, integration_status = NA_character_,
    adaptive_completed = NA, reached_cap = NA_integer_,
    budget_exhausted = NA_integer_, fallback_families = NA_integer_,
    hard_failures = NA_integer_, nonfinite_integrals = NA_integer_,
    max_evaluations = NA_integer_, mode_failures = NA_integer_,
    hessian_regularizations = NA_integer_,
    runtime_seconds = unname(elapsed[["elapsed"]]), log_likelihood = NA_real_,
    fixed_effects = cox_fixed,
    warning_count = as.integer(!is.null(cox_error)), n = nrow(family), families = length(ids),
    max_family_size = max(table(family$famID)), stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE)
print(result)
