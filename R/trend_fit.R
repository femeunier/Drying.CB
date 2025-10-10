trend_fit <- function(v, tn, t0) {
  if (all(!is.finite(v))) return(c(NA_real_, NA_real_))
  ok <- is.finite(v) & is.finite(tn)
  if (sum(ok) < 2) return(c(NA_real_, NA_real_))
  X <- cbind(1, tn[ok] - t0)
  beta <- tryCatch(qr.solve(X, v[ok]), error = function(e) c(NA_real_, NA_real_))
  if (any(!is.finite(beta))) return(c(NA_real_, NA_real_))
  beta  # c(intercept_at_t0, slope_per_year)
}
