convertPressure <- function(barpress, unit){
  # Spit back atm
  if (unit == "atm") {
   barpress_atm <- barpress
  }
  # hPa to atm
  if (unit == "hPa") {
   barpress_atm <- barpress * 0.00098692316931427
  }
  # Torr to atm
  if (unit == "Torr" | unit == "mmHg") {
   barpress_atm <- barpress / 760
  }
  # inHg to atm
  if (unit == "inHg"){
    barpress_atm <- barpress / 29.921
  }
  # psi to atm
  if (unit == "psi") {
   barpress_atm <- barpress * 14.6959487755142
  }
  # bar to atm
  if (unit == "bar") {
   barpress_atm <- barpress * 1.01325
  }
  # kPa to atm
  if(unit == "kPa") {
    barpress_atm <- barpress / 101.325
  }
  # stop message for non-sanctioned unit
  if (!(unit %in% c("atm", "hPa", "Torr", "psi", "bar", "mmHg", "inHg", "kPa"))) {
    stop("Please report barometric pressure in units of `atm`, `hPa`, `psi`,
            `bar`, `mmHg`, `kPa`, or `Torr`.")
  }
  return(barpress_atm)
}
