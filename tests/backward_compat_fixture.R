args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  message("Cross-version fixture skipped: no output path supplied.")
  quit(save = "no", status = 0L)
}
stopifnot(length(args) == 1L)

suppressPackageStartupMessages({
  library(survival)
  library(frailtypack)
})

fixture <- data.frame(
  t0 = 0,
  time = c(1, 2, 3, 4, 1.2, 2.2, 3.2, 4.2, 1.4, 2.4, 3.4, 4.4),
  status = rep(c(1, 0, 1, 0), 3),
  x = rep(c(-0.5, 0.2, 0.7, -0.1), 3),
  fam = rep(1:3, each = 4)
)
K <- kronecker(diag(3), matrix(0.25, 4, 4))
diag(K) <- 1

fit <- suppressWarnings(frailtyPenal(
  Surv(t0, time, status) ~ x + cluster(fam),
  data = fixture,
  hazard = "Weibull",
  RandDist = "LogN",
  recurrentAG = TRUE,
  covMatrix1 = K,
  print.times = FALSE,
  maxit = 1L,
  init.B = 0
))

saveRDS(list(
  version = as.character(packageVersion("frailtypack")),
  b = fit$b,
  logLik = fit$logLik
), args[[1L]])
