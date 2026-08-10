helper <- getFromNamespace(".prepare_two_covariance_init", "frailtypack")

stopifnot(
  isTRUE(all.equal(helper(0.49, 0.25), 0.7)),
  isTRUE(all.equal(helper(c(0.25, 0.49), 0.25), 0.7)),
  isTRUE(all.equal(helper(0, 0), 0))
)

expect_error <- function(expr, pattern) {
  value <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = identity
  )
  stopifnot(inherits(value, "error"), grepl(pattern, conditionMessage(value)))
}

expect_error(helper(c(0.3, 0.49), 0.25), "must equal sig.fixed")
expect_error(helper(-0.1, 0.25), "finite and non-negative")
expect_error(helper(0.49, NULL), "sig.fixed")
