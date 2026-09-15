#' Calculate dissolved oxygen concentrations at saturation given water temperature, atmospheric pressure, and salinity
#'
#' @param temp Water temperature in degrees C
#' @param pressure Barometric pressure, where `pressUnits` denotes units
#' @param pressUnits Units of barometric pressure, must be "atm", "hPa", "Torr", "psi", "bar", "inHg", or "mmHg"
#' @param salinity Water salinity in units of per mille. Default is 0.
#' @param outUnits Preferred units for output, must be "mg/L" or "umol/L"
#'
#' @return A numeric vector
#'
#' @references Garcia, H. E., & Gordon, L. I. (1992). Oxygen solubility in seawater: Better fitting equations. Limnology and Oceanography, 37(6), 1307–1312. https://doi.org/10.4319/lo.1992.37.6.1307
#'
#' @export
#'
#' @examples
#' o2sat(temp = 25, pressure = 29.8, pressUnits = "inHg", outUnits = "mg/L")
o2sat <- function(temp, pressure, pressUnits, salinity = 0, outUnits){
  # Convert barometric pressure into units of atm
  barpress_atm <- convertPressure(barpress = pressure, unit = pressUnits)

  # Vapor pressure correction
  pressCorr <- pressureCorrection(temp, barpress_atm)

  # O2 saturation calculation Combined fit coefficients [umol/kg]
  # (Garcia and Gordon 1992, Table 1)
  A0 <- 5.80818
  A1 <- 3.20684
  A2 <- 4.1189
  A3 <- 4.93845
  A4 <- 1.01567
  A5 <- 1.41575
  B0 <- -7.01211 * 10^-3
  B1 <- -7.25958 * 10^-3
  B2 <- -7.93334 * 10^-3
  B3 <- -5.54491 * 10^-3
  C0 <- -1.32412 * 10^-7
  # Scaled temperature (Garcia and Gordon 1992, eqn. 8)
  TS <- log((298.15 - temp)/(273.15 + temp))  # log() == natural log (ln)
  # Salinity [per mille]
  S <- salinity

  # Calculate O2 saturation concentration at temperature and salinity
  # (Garcia and Gordon 1992, eqn. 8)
  lnO2sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^2 + A3 * TS^3 + A4 *
    TS^4 + A5 * TS^5 + S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3) + C0 * S^2
  O2sat <- exp(lnO2sat)

  # Correct O2 saturation with pressure correction, solubility.conc units
  # [umol/kg]
  result_umolkg <- O2sat * pressCorr

  # Converting from umol/kg to umol/L requires knowing water density at temp
  dens_kgL <- h2odensity(temp)
  result_umolL <- result_umolkg * dens_kgL

  # Convert from umol/L to mg/L
  result_mgL <- result_umolL * 32 / 1000

  # Return desired output to user
  if(outUnits == "umol/L"){
    return(result_umolL)
  }
  if(outUnits == "mg/L"){
    return(result_mgL)
  }
}
