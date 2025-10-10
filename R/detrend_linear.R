detrend_linear <- function(v, tn) {
  if (all(!is.finite(v))) return(v)
  ok <- is.finite(v) & is.finite(tn)
  if (sum(ok) < 2) return(v)
  t0 <- mean(tn[ok])
  # fit v ~ 1 + (tn - t0)
  X <- cbind(1, tn[ok] - t0)
  beta <- tryCatch(qr.solve(X, v[ok]), error = function(e) c(NA, NA))
  if (any(!is.finite(beta))) return(v)
  # remove only the slope component to preserve mean level
  v[ok] <- v[ok] - beta[2] * (tn[ok] - t0)
  v
}
