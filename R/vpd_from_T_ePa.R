vpd_from_T_ePa <- function(T_C, e_Pa) {
  # e_Pa: actual vapour pressure [Pa]
  es_kPa <- 0.6108 * exp((17.27 * T_C) / (T_C + 237.3))
  e_kPa  <- e_Pa / 1000

  pmax(es_kPa - e_kPa, 0)
}
