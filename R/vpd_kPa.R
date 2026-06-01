vpd_kPa <- function(T, Td) {
  es <- function(x)
    0.6108 * exp(17.27 * x / (x + 237.3))
  pmax(0, es(T) - es(Td))
}
