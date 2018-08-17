#' \code{mimsy} functionshortdescription
#'
#' function long description
#'
#' @param filename the name of the .csv file containing raw MIMS data. See github readme for .csv formatting guide
#' @param bgcorr a logical value indicating whether a background correction should be performed. Defaults to FALSE. If TRUE, user must specify values to use for the background correction.
#' @param barpress the barometric pressure of the room while samples were run on the MIMS. User may input a vector, script will take mean value.
#' @param barpress_units units of barometric pressure. Must be one of atm, hPa, psi, bar, or Torr
#' @param salinity salinity of standards, in units of per mille. defaults to 0.
#' @param std.temp a numberic vector (maximum length 2) of the water temperatures (degC) of the standard baths.
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


library(dplyr)

filename <- "MIMS_KAWN24HR_25Mar18_8thSt.csv"
bgcorr <- FALSE
barpress <- 980.35
barpress_units <- "hPa"
std.temp <- c(9.81, 12.05)
salinity <- 0

mimsy <- function(filename, bgcorr = FALSE, barpress, barpress_units, salinity = 0,
                  std.temp){
  #
  #### Raw data import #################################################################
  data <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE) # import csv file
  #
  # UPDATEFLAG: check csv is correct format
  #
  #### Format time column ################################################################
  # as.POSIXct will automatically add in a date to the time that's passed to it, generally this won't be an issue but if samples were analyzed spanning two days (e.g. start at 8pm, end at 1am) this will cause an issue. probably a good idea to come up with a solution. insert check to see if end time > start time, and return warning, for now
  # easiest fix may be dateRun column and paste
  data$Time <- as.POSIXct(data$Time, format = "%I:%M:%S %p")
  #
  if (data$Time[1] > tail(data$Time, n=1)){
    stop("Did you run your samples overnight (e.g. start at 8:00:00 PM, end at 1:00:00 AM the next
day?) Currently mimsy isn't set up to handle this.
\nPlease let me know that this is an issue for you and I'll look into a bug patch asap!")
  }
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
  #### Barometric pressure #############################################################
  #
  if (barpress_units == "atm"){
    Barpress.atm <- mean(barpress)
  }
  if (barpress_units == "hPa"){
    Barpress.atm <- mean(barpress)*0.00098692316931427 # conversion from hPa to atm
  }
  if (barpress_units == "Torr"){
    Barpress.atm <- mean(barpress)*760
  }
  if (barpress_units == "psi"){
    Barpress.atm <- mean(barpress)*14.6959487755142
  }
  if (barpress_units == "bar"){
    Barpress.atm <- mean(barpress)*1.01325
  }
  if (!(barpress_units %in% c("atm", "hPa", "Torr", "psi", "bar"))){
    stop("Please report barometric pressure in units of `atm`, `hPa`, `psi`, `bar`, or `Torr`.")
  }
  #
  #### Calculation of dissolved gasses at saturation ###################################
  #
  # initialize vector to store saturation values
  solubility.conc <- data.frame(O2.conc = numeric(length = 2), N2.conc = numeric(length = 2),
                         Ar.conc = numeric(length = 2), row.names = std.temp)
  #
  # O2 saturation calculation
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
    # Correct O2 saturation with pressure correction, solubility.conc units [umol/kg]
    solubility.conc$O2.conc[i] <- O2.sat * press.corr
  }
  #
  # N2 saturation calculation
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
    # Correct saturation with pressure correction, solubility.conc units are [umol/kg]
    solubility.conc$N2.conc[i] <- N2.sat * press.corr
  }
  #
  # Ar saturation calculation
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
    # Correct saturation with pressure correction, solubility.conc units are [umol/kg]
    solubility.conc$Ar.conc[i] <- Ar.sat * press.corr
  }
  #
  #### Calibration factors ################################################################
  #
  # Group data by Type (`Standard` or `Sample`) and Group (numeric, 1:n)
  data.byTypeGroup <- dplyr::group_by(data, Type, Group)
  #
  # two-point calibration has 6 standard readings
  if (all(data.byTypeGroup[1:6,1] ==
      c("Standard","Standard","Standard","Standard","Standard","Standard"))){
    message("Calculating dissolved concentrations based on a two-point temperature standard.")
    #
    # assemble empty data frame for calibration factors
    calfactor <-
    data.frame(calfactor_28 = numeric(length = max(data.byTypeGroup$Group)*2),
               calfactor_32 = numeric(length = max(data.byTypeGroup$Group)*2),
               calfactor_40 = numeric(length = max(data.byTypeGroup$Group)*2),
               row.names = paste0("STD", rep(1:max(data.byTypeGroup$Group), each = 2), "_",
                                  std.temp))
    #
    for (groupNo in 1:max(data.byTypeGroup$Group)){
      # individually extract each group of standards
      cal.block <-
        data.byTypeGroup %>%
        filter(Type == "Standard" && Group == groupNo)
      #
      # calfactor = solubility.concStd1 / avg(dataStd1_A)
      # Mass28 (N2)
      calfactor$calfactor_28[2*groupNo-1] <- solubility.conc$N2.conc[1] / mean(cal.block$X28[1:3])
      calfactor$calfactor_28[2*groupNo] <- solubility.conc$N2.conc[2] / mean(cal.block$X28[4:6])
      # Mass32 (O2)
      calfactor$calfactor_32[2*groupNo-1] <- solubility.conc$O2.conc[1] / mean(cal.block$X32[1:3])
      calfactor$calfactor_32[2*groupNo] <- solubility.conc$O2.conc[2] / mean(cal.block$X32[4:6])
      # Mass40 (Ar)
      calfactor$calfactor_40[2*groupNo-1] <- solubility.conc$Ar.conc[1] / mean(cal.block$X40[1:3])
      calfactor$calfactor_40[2*groupNo] <- solubility.conc$Ar.conc[2] / mean(cal.block$X40[4:6])
    }
    #
    # calculate drift slope
    # add drift slope column to calfactor
    calfactor$driftslope_28 <- NA
    calfactor$driftslope_32 <- NA
    calfactor$driftslope_40 <- NA
    # subset dataframe of just the standard readings
    data.standards <-
      data.byTypeGroup %>%
      filter(Type == "Standard")
    #
    # perform drift slope calculations
    for (groupNo in 1:max(data.byTypeGroup$Group)){
      #
      # set up a vector for the time indexes, temp1
      time.indexer <- numeric(length = max(data.byTypeGroup$Group))
      time.indexer[1] <- 7
      for(i in 2:length(time.indexer)){
        time.indexer[i] <- time.indexer[i-1] + 6
      }
      # corrections for temp 1
      # Mass28
      calfactor$driftslope_28[2*groupNo-1] <-
        (calfactor$calfactor_28[2*groupNo+1] - calfactor$calfactor_28[2*groupNo-1]) /
        as.numeric(difftime(data.standards$Time[time.indexer[groupNo]],
                            data.standards$Time[2*groupNo-1], units = "mins"))
      # Mass32
      calfactor$driftslope_32[2*groupNo-1] <-
        (calfactor$calfactor_32[2*groupNo+1] - calfactor$calfactor_32[2*groupNo-1]) /
        as.numeric(difftime(data.standards$Time[time.indexer[groupNo]],
                            data.standards$Time[2*groupNo-1], units = "mins"))
      # Mass40
      calfactor$driftslope_40[2*groupNo-1] <-
        (calfactor$calfactor_40[2*groupNo+1] - calfactor$calfactor_40[2*groupNo-1]) /
        as.numeric(difftime(data.standards$Time[time.indexer[groupNo]],
                            data.standards$Time[2*groupNo-1], units = "mins"))
      #
      # set up a vector for the time indicies, temp2
      time.indexer[1] <- 10
      for(i in 2:length(time.indexer)){
        time.indexer[i] <- time.indexer[i-1] + 6
      }
      #
      # corrections for temp 2
      # Mass28
      calfactor$driftslope_28[2*groupNo] <-
        (calfactor$calfactor_28[2*groupNo+2] - calfactor$calfactor_28[2*groupNo]) /
        as.numeric(difftime(data.standards$Time[time.indexer[groupNo]],
                            data.standards$Time[2*groupNo+2], units = "mins"))
      # Mass32
      calfactor$driftslope_32[2*groupNo] <-
        (calfactor$calfactor_32[2*groupNo+2] - calfactor$calfactor_32[2*groupNo]) /
        as.numeric(difftime(data.standards$Time[time.indexer[groupNo]],
                            data.standards$Time[2*groupNo+2], units = "mins"))
      # Mass40
      calfactor$driftslope_40[2*groupNo] <-
        (calfactor$calfactor_40[2*groupNo+2] - calfactor$calfactor_40[2*groupNo]) /
        as.numeric(difftime(data.standards$Time[time.indexer[groupNo]],
                          data.standards$Time[2*groupNo+2], units = "mins"))
    }
    #
    # drift corrected cal factor: calibration factor + (slope * (standard start time - time[i]))
    #
    # add columns to data.byTypeGroup
    data.byTypeGroup$CalFactor_28 <- numeric(length = nrow(data.byTypeGroup))
    data.byTypeGroup$CalFactor_32 <- numeric(length = nrow(data.byTypeGroup))
    data.byTypeGroup$CalFactor_40 <- numeric(length = nrow(data.byTypeGroup))
    #
    #
    for(groupNo in 1:max(data.byTypeGroup$Group)){
      # want to deal with one sampling group at a time
      sampleGroup <-
      data.byTypeGroup %>%
        filter(Group == groupNo)
      #
      for(row in 1:length(sampleGroup)){
        # calculate calibration factor for each sample reading, the calibration factor is calculated based only on the initial temperature standard (as suggested in MIMS documentation)
        # UPDATEFLAG: check this ^^^
        # temp 1
        # Mass28
        data.byTypeGroup$CalFactor_28[row] <- # change this indexing, row will only go 1-12
        calfactor$calfactor_28[2*groupNo-1] +
          (calfactor$driftslope_28[2*groupNo-1] *
             (as.numeric(difftime(sampleGroup$Time[row], sampleGroup$Time[1], units = "mins"))))
        # Mass32
        data.byTypeGroup$CalFactor_32[row] <-
          calfactor$calfactor_32[2*groupNo-1] +
          (calfactor$driftslope_32[2*groupNo-1] *
             (as.numeric(difftime(sampleGroup$Time[row], sampleGroup$Time[1], units = "mins"))))
        # Mass40
        data.byTypeGroup$CalFactor_28[row] <-
          calfactor$calfactor_28[2*groupNo-1] +
          (calfactor$driftslope_28[2*groupNo-1] *
             (as.numeric(difftime(sampleGroup$Time[row], sampleGroup$Time[1], units = "mins"))))
      }
    }
    # get a concentration value in microM = raw MIMS reading * drift corr cal
    #
    # return list solubility.conc, calfactors, longformdata (data), shortformdata (just relavant units)
  }
  #
  #
  # single point calibration #############################################################
  if (all(data.byTypeGroup[1:6,1] ==
      c("Standard","Standard","Standard","Sample","Sample","Sample"))){
    # add in error message for now
    stop("Single-point temperature calibration not yet supported")
    #
    message("Calculating dissolved concentrations based on a one point temperature standard.
            \nIf this is an error, check if standard readings have been properly designated `Standard` in the `Type` column.")
    #
    calfactor <-
    data.frame(calfactor_28 = numeric(length = max(data.byTypeGroup$Group)),
               calfactor_32 = numeric(length = max(data.byTypeGroup$Group)),
               calfactor_40 = numeric(length = max(data.byTypeGroup$Group)),
               row.names = paste0("STD", 1:max(data.byTypeGroup$Group))
    )
    }
}









