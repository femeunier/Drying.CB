vpd_from_T_q <- function(T_C, q, p_kPa) {
  # T_C: air temperature in °C (vector)
  # q: specific humidity in kg/kg (vector)
  # p_kPa: air pressure in kPa (vector or scalar; ~101.325 at sea level)
  eps <- 0.622
  es_kPa <- 0.6108 * exp((17.27 * T_C) / (T_C + 237.3))
  e_kPa  <- (q * p_kPa) / (eps + (1 - eps) * q)
  pmax(es_kPa - e_kPa, 0)
}
