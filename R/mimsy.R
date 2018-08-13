#' \code{mimsy} functionshortdescription
#'
#' function long description
#'
#' @param filename the name of the .csv file containing raw MIMS data. See github readme for .csv formatting guide
#' @param bgcorr a logical value indicating whether a background correction should be performed. Defaults to FALSE. If TRUE, user must specify values to use for the background correction.
#' @param barpress_units units of measured atmospheric pressure. function recognizes "atm", "hPa", (more units coming with update...)
#' @param salinity salinity of standards, in units of per mille. defaults to 0.
#' @param twoPoint a logical value indicating whether a two point temperature calibration was used (TRUE) or a single point temperature calibration was used (FALSE). Defaults to TRUE.
#' @param std.temp a numberic vector (maximum length 2) of the water temperatures used for the standard calibration.
#'
#' @return dataframe of corrected sample readings in units of microM or mg
#'
#' @author Michelle Catherine Kelly
#'
#' @references
#'
#' Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater: Better fitting
#' equations}, Limnology and Oceanography, 37(6).
#'
#' Benson, B. B. & Krause, D. (1984). \emph{The concentration and isotopic
#' fractionation of oxygen dissolved in freshwater and seawater in equilibrium
#' with the atmosphere.} Limnology and Oceanography, 29(3), 620-632.
#' doi:10.4319/lo.1984.29.3.0620
#'
#' Stull, D. R. (1947). \emph{Vapor Pressure of Pure Substances. Organic and
#' Inorganic Compounds.} Industrial & Engineering Chemistry, 39(4), 517-540.
#' doi: 10.1021/ie50448a022
#'
#' Hamme, R. C. & Emerson, S. R. (2004). \emph{The solubility of neon, nitrogen and argon
#' in distilled water and seawater}, Deep-Sea Research I, 51(11), 1517-1528.
#'
#' @examples
#' mydata <- mimsy(filename = "./RawData/MIMS_data.csv", bgcorr = FALSE, std.temp = c(12.5, 13.2))
#'
#' @export

#filename <- "MIMS_KAWN24HR_25Mar18_8thSt.csv"
#bgcorr <- FALSE
#barpress_units <- "hPa"
#std.temp <- c(9.81, 12.05)
#salinity <- 0

mimsy <- function(filename, bgcorr = FALSE, barpress_units, salinity = 0,
                  twoPoint = TRUE, std.temp){
  #
  #### Raw data import #################################################################
  data <- read.csv(filename, header = TRUE) # import csv file
  #
  # UPDATEFLAG: check csv is correct format
  # UPDATEFLAG: check that run date is the same, csv only contains one run date of data
  #
  #### Barometric pressure #############################################################
  #
  # UPDATEFLAG: add more unit conversions. check with weather stations on their units
  if (barpress_units == "atm"){
    Barpress.atm <- mean(data$Barpress)
  }
  if (barpress_units == "hPa"){
    Barpress.atm <- mean(data$Barpress)*0.00098692316931427 # conversion from hPa to atm
  }
  #
  #
  #
  #### Calculation of dissolved gasses at saturation ###################################
  #
  # initialize vector to store saturation values
  sat.conc <- data.frame(O2.conc = numeric(length = 2), N2.conc = numeric(length = 2),
                         Ar.conc = numeric(length = 2), row.names = std.temp)
  #
  # O2 saturation calculation ##########################################################
  for (i in seq_along(std.temp)) {
    t <- std.temp[i]
    #### Vapor pressure correction ####
    # use the Antoine equation to calculate vapor pressure of water [bar]
    # See NIST Chemistry WebBook for general tables, these parameters valid for
    # temperatures between -18 to 100C (Stull 1947)
    vapor.press <- exp(4.6543 - (1435.264 / ((t+273.15) + -64.848)))
    vapor.press <- vapor.press * 0.98692 # conversion from [bar] to [atm]
    #
    # pressure correction [atm] =
    #       (current pressure - vapor pressure) / (standard pressure [atm] - vapor pressure)
    press.corr <- (Barpress.atm - vapor.press) / (1 - vapor.press)
    #
    #### O2 saturation calculation ####
    # Combined fit coefficients [umol/kg] (Garcia and Gordon 1992, Table 1)
    A0 <- 5.80818
    A1 <- 3.20684
    A2 <- 4.11890
    A3 <- 4.93845
    A4 <- 1.01567
    A5 <- 1.41575
    B0 <- -7.01211*10^-3
    B1 <- -7.25958*10^-3
    B2 <- -7.93334*10^-3
    B3 <- -5.54491*10^-3
    C0 <- -1.32412*10^-7
    #
    # Scaled temperature (Garcia and Gordon 1992, eqn. 8)
    TS <- log((298.15 - t)/(273.15 + t)) # log() == natural log (ln)
    #
    # Salinity [per mille]
    S <- salinity
    # Calculate O2 saturation concentration at temperature and salinity (Garcia and Gordon
    # 1992, eqn. 8)
    lnO2.sat <- A0 + A1*TS + A2*TS^2 + A3*TS^2 +
      A3*TS^3 + A4*TS^4 + A5*TS^5 +
      S*(B0 + B1*TS + B2*TS^2 + B3*TS^3) +
      C0*S^2
    O2.sat <- exp(lnO2.sat)
    # Correct O2 saturation with pressure correction, sat.conc units [umol/kg]
    sat.conc$O2.conc[i] <- O2.sat * press.corr
  }
  #
  # N2 saturation calculation ##############################################################
  for (i in seq_along(std.temp)) {
    t <- std.temp[i]
    #### Vapor pressure correction ####
    # use the Antoine equation to calculate vapor pressure of water [bar]
    # See NIST Chemistry WebBook for general tables, these parameters valid for
    # temperatures between -18 to 100C (Stull 1947)
    vapor.press <- exp(4.6543 - (1435.264 / ((t+273.15) + -64.848)))
    vapor.press <- vapor.press * 0.98692 # conversion from [bar] to [atm]
    #
    # pressure correction [atm] =
    #       (current pressure - vapor pressure) / (standard pressure [atm] - vapor pressure)
    press.corr <- (Barpress.atm - vapor.press) / (1 - vapor.press)
    #
    #### N2 saturation calculation ####
    # Coefficients [umol/kg] (Hamme and Emerson 2004, Table 4)
    A0 <- 6.42931
    A1 <- 2.92704
    A2 <- 4.32531
    A3 <- 4.69149
    B0 <- -7.44129*10^-3
    B1 <- -8.02566*10^-3
    B2 <- -1.46775*10^-2
    # Scaled temperature (Hamme and Emerson 2004, eqn. 2, but identical to
    # Garcia and Gordon 1992, eqn. 8)
    TS <- log((298.15 - t)/(273.15 + t))
    #
    # Salinity [per mille]
    S <- salinity
    # Calculate saturation concentration at temperature and salinity (Hamme and Emerson 2004,
    # eqn. 1)
    lnN2.sat <- A0 + A1*TS + A2*TS^2 + A3*TS^3 +
      S*(B0 + B1*TS + B2*TS^2 + B3*TS^3)
    N2.sat <- exp(lnN2.sat)
    # Correct saturation with pressure correction, sat.conc units are [umol/kg]
    sat.conc$N2.conc[i] <- N2.sat * press.corr
  }
  #
  # Ar saturation calculation ##############################################################
  for (i in seq_along(std.temp)) {
    t <- std.temp[i]
    #### Vapor pressure correction ####
    # use the Antoine equation to calculate vapor pressure of water [bar]
    # See NIST Chemistry WebBook for general tables, these parameters valid for
    # temperatures between -18 to 100C (Stull 1947)
    vapor.press <- exp(4.6543 - (1435.264 / ((t+273.15) + -64.848)))
    vapor.press <- vapor.press * 0.98692 # conversion from [bar] to [atm]
    #
    # pressure correction [atm] =
    #       (current pressure - vapor pressure) / (standard pressure [atm] - vapor pressure)
    press.corr <- (Barpress.atm - vapor.press) / (1 - vapor.press)
    #
    #### Ar saturation calculation ####
    # Coefficients [umol/kg] (Hamme and Emerson 2004, Table 4)
    A0 <- 2.79150
    A1 <- 3.17609
    A2 <- 4.13116
    A3 <- 4.90379
    B0 <- -6.96233*10^-3
    B1 <- -7.99970*10^-3
    B2 <- -1.16888*10^-2
    # Scaled temperature (Hamme and Emerson 2004, eqn. 2, but identical to
    # Garcia and Gordon 1992, eqn. 8)
    TS <- log((298.15 - t)/(273.15 + t))
    #
    # Salinity [per mille]
    S <- salinity
    # Calculate saturation concentration at temperature and salinity (Hamme and Emerson 2004,
    # eqn. 1)
    lnAr.sat <- A0 + A1*TS + A2*TS^2 + A3*TS^3 +
      S*(B0 + B1*TS + B2*TS^2 + B3*TS^3)
    Ar.sat <- exp(lnAr.sat)
    # Correct saturation with pressure correction, sat.conc units are [umol/kg]
    sat.conc$Ar.conc[i] <- Ar.sat * press.corr
  }
  #
  #
  #### Background corrections ###############################################################
  #
  if (bgcorr == FALSE) {
    bgcorr <- 0 # set background correction value to zero
  }
  if (bgcorr != FALSE){ # UPDATEFLAG
    message("Background correction not yet supported. In meantime, please set bgcorr to FALSE")
  }
  #
  #### Calibration factors
  dplyr::group_by
}









