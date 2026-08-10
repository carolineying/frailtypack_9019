#!/usr/bin/env Rscript
files <- commandArgs(trailingOnly = TRUE)
if (!length(files)) stop("Supply one or more TSV files from compare_adaptive_genz_keister.R")
x <- do.call(rbind, lapply(files, read.delim, check.names = FALSE))
key <- paste(x$method, ifelse(is.na(x$nominal_maxpts), "NA", x$nominal_maxpts),
             sep = ":")
summary <- do.call(rbind, lapply(split(x, key), function(d) data.frame(
  method = d$method[[1L]], maxpts = d$nominal_maxpts[[1L]],
  replicates = nrow(d), mean_sigma2 = mean(d$sigma2_hat, na.rm = TRUE),
  bias = mean(d$bias, na.rm = TRUE),
  rmse = sqrt(mean(d$squared_error, na.rm = TRUE)),
  convergence_rate = mean(d$converged, na.rm = TRUE),
  valid_result_rate = mean(d$result_valid, na.rm = TRUE),
  integration_failure_rate = mean(d$integration_failure, na.rm = TRUE),
  cap_rate = mean(d$reached_cap > 0, na.rm = TRUE),
  budget_exhaustion_rate = mean(d$budget_exhausted > 0, na.rm = TRUE),
  fallback_rate = mean(d$fallback_families > 0, na.rm = TRUE),
  hard_failure_rate = mean(d$hard_failures > 0, na.rm = TRUE),
  nonfinite_integral_rate = mean(d$nonfinite_integrals > 0, na.rm = TRUE),
  mean_runtime_seconds = mean(d$runtime_seconds, na.rm = TRUE),
  ordering_max_abs_difference = {
    ranges <- vapply(split(d$sigma2_hat, d$replicate), function(v)
      if (length(v) > 1L) diff(range(v, na.rm = TRUE)) else NA_real_, numeric(1))
    if (all(is.na(ranges))) NA_real_ else max(ranges, na.rm = TRUE)
  },
  mode_failures = sum(d$mode_failures, na.rm = TRUE),
  hessian_regularizations = sum(d$hessian_regularizations, na.rm = TRUE)
)))
print(summary, row.names = FALSE)
