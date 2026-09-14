pressureCorrection <- function(temp, barpress_atm){
  # Use the Antoine equation to calculate vapor pressure of water [bar] See
  # NIST Chemistry WebBook for general tables, these parameters valid for
  # temperatures between -18 to 100C (Stull 1947)
  vaporPress <- exp(4.6543 - (1435.264/((temp + 273.15) + -64.848)))
  vaporPress <- vaporPress * 0.98692  # conversion from [bar] to [atm]

  # pressure correction [atm] = (current pressure - vapor pressure) /
  # (standard pressure [atm] - vapor pressure)
  pressCorr <- (barpress_atm - vaporPress)/(1 - vaporPress)
  return(pressCorr)
}
