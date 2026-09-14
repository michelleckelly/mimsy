#' Calculate argon gas concentrations at saturation given water temperature, atmospheric pressure, and salinity
#'
#' @param temp Water temperature in degrees C
#' @param pressure Barometric pressure, where `pressUnits` denotes units
#' @param pressUnits Units of barometric pressure, must be "atm", "hPa", "Torr", "psi", "bar", "inHg", or "mmHg"
#' @param salinity Water salinity in units of per mille. Default is 0.
#' @param outUnits Preferred units for output, must be "mg/L" or "umol/L"
#'
#' @return A numeric vector
#'
#' @references Hamme, R. C., & Emerson, S. R. (2004). The solubility of neon, nitrogen and argon in distilled water and seawater. Deep Sea Research Part I: Oceanographic Research Papers, 51(11), 1517–1528. https://doi.org/10.1016/j.dsr.2004.06.009
#'
#' @export
#'
#' @examples
#' arsat(temp = 25, pressure = 29.8, pressUnits = "inHg", outUnits = "mg/L")
arsat <- function(temp, pressure, pressUnits, salinity = 0, outUnits){
  # Convert barometric pressure into units of atm
  barpress_atm <- convertPressure(barpress = pressure, unit = pressUnits)

  # Vapor pressure correction
  pressCorr <- pressureCorrection(temp, barpress_atm)

  # Ar saturation calculation Coefficients [umol/kg]
  # (Hamme and Emerson 2004, Table 4)
  A0 <- 2.7915
  A1 <- 3.17609
  A2 <- 4.13116
  A3 <- 4.90379
  B0 <- -6.96233 * 10^-3
  B1 <- -7.9997 * 10^-3
  B2 <- -1.16888 * 10^-2
  # Scaled temperature (Hamme and Emerson 2004, eqn. 2,
  # but identical to Garcia and Gordon 1992, eqn. 8)
  TS <- log((298.15 - temp)/(273.15 + temp))

  # Salinity [per mille]
  S <- salinity

  # Calculate saturation concentration at temperature and salinity
  # (Hamme and emerson 2004, eqn. 1)
  lnArsat <-
    A0 + A1 * TS + A2 * TS^2 + A3 * TS^3 + S * (B0 + B1 * TS + B2 * TS^2)
  Arsat <- exp(lnArsat)

  # Correct N2 saturation with pressure correction, solubility.conc units
  # [umol/kg]
  result_umolkg <- Arsat * pressCorr

  # Converting from umol/kg to umol/L requires knowing water density at temp
  dens_kgL <- h2odensity(temp)
  result_umolL <- result_umolkg * dens_kgL

  # Convert from umol/L to mg/L
  result_mgL <- result_umolL * 40 / 1000

  # Return desired output to user
  if(outUnits == "umol/L"){
    return(result_umolL)
  }
  if(outUnits == "mg/L"){
    return(result_mgL)
  }
}
