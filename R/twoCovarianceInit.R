.prepare_two_covariance_init <- function(init.Theta, sig.fixed) {
  if (!is.numeric(sig.fixed) || length(sig.fixed) != 1L ||
      !is.finite(sig.fixed) || sig.fixed < 0) {
    stop("sig.fixed must be one finite, non-negative variance when two covariance matrices are provided.")
  }

  if (!is.numeric(init.Theta) || !length(init.Theta) %in% c(1L, 2L) ||
      any(!is.finite(init.Theta)) || any(init.Theta < 0)) {
    stop(
      paste(
        "For two covariance matrices, init.Theta must contain either the",
        "free covMatrix2 variance or c(sig.fixed, free variance), with all",
        "values finite and non-negative."
      )
    )
  }

  if (length(init.Theta) == 2L &&
      !isTRUE(all.equal(as.numeric(init.Theta[[1L]]), as.numeric(sig.fixed),
                        tolerance = sqrt(.Machine$double.eps)))) {
    stop("The first element of a two-value init.Theta must equal sig.fixed.")
  }

  free_variance <- as.numeric(init.Theta[[length(init.Theta)]])

  # The compiled likelihood represents a variance as an unconstrained
  # parameter squared. Keep the public interface on the documented variance
  # scale and convert only at the R/Fortran boundary.
  sqrt(free_variance)
}
