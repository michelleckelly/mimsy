#' Calculate water density at temperatures between 0 and 40 degC
#'
#' @param temp Water temperature in degrees C
#'
#' @return A numeric vector in units of kg/m3
#'
#' @references Tanaka, M., Girard, G., Davis, R., Peuto, A., & Bignell, N. (2001). Recommended table for the density of water between 0 C and 40 C based on recent experimental reports. Metrologia, 38(4), 301–309. https://doi.org/10.1088/0026-1394/38/4/3
#'
#'
#' @export
h2odensity <- function(temp){

  # Fit coefficients and equation from Tanaka et al 2001
  a1 <- -3.983035
  a2 <- 301.797
  a3 <- 522528.9
  a4 <- 69.34881
  a5 <- 999.974950 # kg/m3

  # Equation on pg 305, Tanaka et al 2001
  numer <- (temp + a1)^2 * (temp + a2)
  denom <- a3 * (temp + a4)

  rho_kgm3 <- a5 * (1 - numer/denom)

  # Convert from kg/m3 to kg/L
  rho_kgL <- rho_kgm3 / 1000
  return(rho_kgL)
}
