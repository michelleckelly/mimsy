#' Calculate mass concentrations in mg/L given molar concentrations in umol/L
#'
#' @param species Molecule, must be "N2", "O2", or "Ar"
#' @param value Concentration of molecule, numeric value
#' @param inUnits Units of molecule concentration, default is "umol/L"
#' @param outUnits Desired units for output, must be "mg/L", "g/L", or "ug/L". default is "mg/L"
#'
#' @return A numeric vector
#'
#' @export
#'
#' @examples
#' molar_to_mass(species = "N2", value = 12.1, inUnits = "umol/L", outUnits = "mg/L")
molar_to_mass <- function(species, value, inUnits = "umol/L", outUnits = "mg/L"){
  # Convert from incoming units to base molar
  if(inUnits == "umol/L"){
    value_molL <- value * 10^-6
  }

  # Convert from molar mass
  if(species == "N2"){
    value_gL <- value_molL * 28.014
  }
  if(species == "O2"){
    value_gL <- value_molL * 31.9988
  }
  if(species == "Ar"){
    value_gL <- value_molL * 39.948
  }

  # Convert to outgoing units
  if(outUnits == "mg/L"){
    outVal <- value_gL * 10^3
  }
  if(outUnits == "g/L"){
    outVal <- value_gL
  }
  if(outUnits == "ug/L"){
    outVal <- value_gL * 10^6
  }

  # Output
  return(outVal)
}
