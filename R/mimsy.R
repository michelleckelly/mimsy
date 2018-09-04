#' \code{mimsy} functionshortdescription
#'
#' function long description
#'
#' @param filename the name of the .csv file containing raw MIMS data. See github readme for .csv formatting guide
#' @param bgcorr a logical value indicating whether a background correction should be performed. Defaults to FALSE. If TRUE, user must specify values to use for the background correction.
#' @param barpress the barometric pressure of the room while samples were run on the MIMS. User may input a vector, script will take mean value.
#' @param barpress_units units of barometric pressure. Must be one of atm, hPa, psi, bar, or Torr
#' @param salinity salinity of standards, in units of per mille. defaults to 0.
#' @param std.tempss a numberic vector (maximum length 2) of the water temperatures (degC) of the standard baths.
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
#' mydata <- mimsy(filename = "./RawData/MIMS_data.csv", bgcorr = FALSE, std.temps = c(12.5, 13.2))
#'
#' @export
#'

library(magrittr)
library(dplyr)

filename <- "MIMS_KAWN24HR_25Mar18_8thSt_mimsyformat.csv"
barpress <- c(981.2, 979.5)
barpress_units <- "hPa"
std.tempss <- c(9.81, 12.05)

mimsy <- function(filename, bgcorr = FALSE, barpress,
                  barpress_units, salinity = 0, std.tempsss){

  # 1. Raw data import ---------------------------------------------
  data <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)

  # UPDATEFLAG: check csv is correct format

  # Format time column ----------------------------------------------

  # as.POSIXct will automatically add in a date to the time that's passed to it, generally this won't be an issue but if samples were analyzed spanning two days (e.g. start at 8pm, end at 1am) this will cause an issue. probably a good idea to come up with a solution. insert check to see if end time > start time, and return warning, for now
  # easiest fix may be dateRun column and paste
  data$Time <- as.POSIXct(data$Time, format = "%I:%M:%S %p")

  if (data$Time[1] > tail(data$Time, n=1)){
    stop("Did you run your samples overnight (e.g. start at 8:00:00 PM, end at 1:00:00 AM the next
day?) Currently mimsy isn't set up to handle this.
\nPlease let me know that this is an issue for you and I'll look into a bug patch asap!")
  }

  # 2. Background corrections -----------------------------------------------

  if (bgcorr == FALSE) {
    bgcorr <- 0 # set background correction value to zero
  }
  if (bgcorr != FALSE){ # UPDATEFLAG
    message("Background correction not yet supported.
            In meantime, please set bgcorr to FALSE")
  }

  # Barometric pressure conversion ------------------------------------------

  if (barpress_units == "atm"){
    Barpress.atm <- mean(barpress)
  }
  # hPa to atm
  if (barpress_units == "hPa"){
    Barpress.atm <- mean(barpress) * 0.00098692316931427
  }
  # Torr to atm
  if (barpress_units == "Torr"){
    Barpress.atm <- mean(barpress) * 760
  }
  # psi to atm
  if (barpress_units == "psi"){
    Barpress.atm <- mean(barpress) * 14.6959487755142
  }
  # bar to atm
  if (barpress_units == "bar"){
    Barpress.atm <- mean(barpress) * 1.01325
  }
  # stop message for non-sanctioned units
  if (!(barpress_units %in% c("atm", "hPa", "Torr", "psi", "bar"))){
    stop("Please report barometric pressure in units of `atm`, `hPa`, `psi`, `bar`, or `Torr`.")
  }

  # 3. Calculate solubilites of dissolved gas ---------------------------

  # initialize vector to store concentration values
  solubility.conc <- data.frame(O2.conc = numeric(length = 2),
                                N2.conc = numeric(length = 2),
                                Ar.conc = numeric(length = 2), row.names = std.temps)

  # O2 saturation calculation -----------------------------------------------------
  for (i in seq_along(std.temps)) {

    t <- std.temps[i]

    # Vapor pressure correction
    # use the Antoine equation to calculate vapor pressure of water [bar]
    # See NIST Chemistry WebBook for general tables, these parameters valid for
    # temperatures between -18 to 100C (Stull 1947)
    vapor.press <- exp(4.6543 - (1435.264 / ((t + 273.15) + -64.848)))
    vapor.press <- vapor.press * 0.98692 # conversion from [bar] to [atm]

    # pressure correction [atm] =
    #       (current pressure - vapor pressure) / (standard pressure [atm] - vapor pressure)
    press.corr <- (Barpress.atm - vapor.press) / (1 - vapor.press)

    # O2 saturation calculation
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
    lnO2.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^2 +
      A3 * TS^3 + A4 * TS^4 + A5 * TS^5 +
      S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3) +
      C0 * S^2
    O2.sat <- exp(lnO2.sat)

    # Correct O2 saturation with pressure correction, solubility.conc units [umol/kg]
    solubility.conc$O2.conc[i] <- O2.sat * press.corr
    }

  # N2 saturation calculation -------------------------------------------------------
  for (i in seq_along(std.temps)) {

    t <- std.temps[i]

    # Vapor pressure correction
    # use the Antoine equation to calculate vapor pressure of water [bar]
    # See NIST Chemistry WebBook for general tables, these parameters valid for
    # temperatures between -18 to 100C (Stull 1947)
    vapor.press <- exp(4.6543 - (1435.264 / ((t + 273.15) + -64.848)))
    vapor.press <- vapor.press * 0.98692 # conversion from [bar] to [atm]


    # pressure correction [atm] =
    #       (current pressure - vapor pressure) / (standard pressure [atm] - vapor pressure)
    press.corr <- (Barpress.atm - vapor.press) / (1 - vapor.press)

    # N2 saturation calculation
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

    # Salinity [per mille]
    S <- salinity

    # Calculate saturation concentration at temperature and salinity
    # (Hamme and Emerson 2004, eqn. 1)
    lnN2.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^3 +
      S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3)
    N2.sat <- exp(lnN2.sat)

    # Correct saturation with pressure correction, solubility.conc units are [umol/kg]
    solubility.conc$N2.conc[i] <- N2.sat * press.corr
  }

  # Ar saturation calculation --------------------------------------------------------
  for (i in seq_along(std.temps)) {

    t <- std.temps[i]

    # Vapor pressure correction
    # use the Antoine equation to calculate vapor pressure of water [bar]
    # See NIST Chemistry WebBook for general tables, these parameters valid for
    # temperatures between -18 to 100C (Stull 1947)
    vapor.press <- exp(4.6543 - (1435.264 / ((t + 273.15) + -64.848)))
    vapor.press <- vapor.press * 0.98692 # conversion from [bar] to [atm]

    # pressure correction [atm] =
    #       (current pressure - vapor pressure) / (standard pressure [atm] - vapor pressure)
    press.corr <- (Barpress.atm - vapor.press) / (1 - vapor.press)

    # Ar saturation calculation
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

    # Salinity [per mille]
    S <- salinity

    # Calculate saturation concentration at temperature and salinity
    # (Hamme and Emerson 2004, eqn. 1)
    lnAr.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^3 +
      S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3)
    Ar.sat <- exp(lnAr.sat)
    # Correct saturation with pressure correction, solubility.conc units are [umol/kg]
    solubility.conc$Ar.conc[i] <- Ar.sat * press.corr
  }

  # 4. Calculate calibration factors ----------------------------------

  # Group data by Type (`Standard` or `Sample`) and Group (numeric, 1:n)
  data <- dplyr::group_by(data, Type, Group)

  # two-point calibration has 6 standard readings
  if (all(data[1:6,1] ==
      c("Standard","Standard","Standard",
        "Standard","Standard","Standard"))){

    message("Calculating dissolved concentrations based on a two-point temperature standard.")

    # assemble empty data frame for calibration factors
    calfactor <-
    data.frame(calfactor_28 = numeric(length = max(data$Group)*2),
               calfactor_32 = numeric(length = max(data$Group)*2),
               calfactor_40 = numeric(length = max(data$Group)*2),
               row.names = paste0("Group", rep(1:max(data$Group), each = 2), "_",
                                  std.temps))

    for (groupNo in 1:max(data$Group)){
      # individually extract each group of standards
      cal.block <-
        data %>%
        filter(Type == "Standard" && Group == groupNo)

      # calculate calibration factor
      # calfactor = solubility.concStd1 / avg(dataStd1_A)

      # Mass28 (N2)
      calfactor$calfactor_28[2*groupNo-1] <-
        solubility.conc$N2.conc[1] / mean(cal.block$X28[1:3])
      calfactor$calfactor_28[2*groupNo] <-
        solubility.conc$N2.conc[2] / mean(cal.block$X28[4:6])
      # Mass32 (O2)
      calfactor$calfactor_32[2*groupNo-1] <-
        solubility.conc$O2.conc[1] / mean(cal.block$X32[1:3])
      calfactor$calfactor_32[2*groupNo] <-
        solubility.conc$O2.conc[2] / mean(cal.block$X32[4:6])
      # Mass40 (Ar)
      calfactor$calfactor_40[2*groupNo-1] <-
        solubility.conc$Ar.conc[1] / mean(cal.block$X40[1:3])
      calfactor$calfactor_40[2*groupNo] <-
        solubility.conc$Ar.conc[2] / mean(cal.block$X40[4:6])
    }

    # 5. Calculate slope and intercepts of calibration curve -----------------
    calslope <-
      data.frame(calslope_28 = numeric(length = max(data$Group)),
                 calintercept_28 = numeric(length = max(data$Group)),
                 calslope_32 = numeric(length = max(data$Group)),
                 calintercept_32 = numeric(length = max(data$Group)),
                 calslope_40 = numeric(length = max(data$Group)),
                 calintercept_40 = numeric(length = max(data$Group)),
                 row.names = paste0("Group", 1:max(data$Group)))

    for (groupNo in 1:max(data$Group)){

      # using a linear model to find slope and intercept
      # of calibration line
      # linear model: y ~ x

      # Mass28
      lm <- lm( c(calfactor$calfactor_28[2*groupNo-1],
                  calfactor$calfactor_28[2*groupNo]) ~ std.temps)
      calslope$calslope_28[groupNo] <- lm$coefficients[2] # slope
      calslope$calintercept_28[groupNo] <- lm$coefficients[1] #intercept

      # Mass32
      lm <- lm( c(calfactor$calfactor_32[2*groupNo-1],
                  calfactor$calfactor_32[2*groupNo]) ~ std.temps)
      calslope$calslope_32[groupNo] <- lm$coefficients[2] # slope
      calslope$calintercept_32[groupNo] <- lm$coefficients[1] #intercept

      # Mass40
      lm <- lm( c(calfactor$calfactor_40[2*groupNo-1],
                  calfactor$calfactor_40[2*groupNo]) ~ std.temps)
      calslope$calslope_40[groupNo] <- lm$coefficients[2] # slope
      calslope$calintercept_40[groupNo] <- lm$coefficients[1] #intercept
    }

    # 6. Perform drift correction for calibration slope and intercept --------

    # add columns to calslope
    calslope$DRIFT.calslope_28 <- NA
    calslope$DRIFT.calintercept_28 <- NA
    calslope$DRIFT.calslope_32 <- NA
    calslope$DRIFT.calintercept_32 <- NA
    calslope$DRIFT.calslope_40 <- NA
    calslope$DRIFT.calintercept_40 <- NA

    for (groupNo in 1:max(data$Group)){
      # individually extract each group
      group.block <-
        data %>%
        filter(Group == groupNo)

      # using a linear model to find drift slope and intercept

      # Drift corrected Mass28
      lm <- lm(calslope$calslope_28[groupNo:(groupNo+1)] ~
                 c(group.block$Time[1], tail(group.block$Time, n = 1)))
      calslope$DRIFT.calslope_28[groupNo] <- lm$coefficients[2] # slope
      calslope$DRIFT.calintercept_28[groupNo] <- lm$coefficients[1] # intercept

      # Drift corrected Mass32
      lm <- lm(calslope$calslope_32[groupNo:(groupNo+1)] ~
                 c(group.block$Time[1], tail(group.block$Time, n = 1)))
      calslope$DRIFT.calslope_32[groupNo] <- lm$coefficients[2] # slope
      calslope$DRIFT.calintercept_32[groupNo] <- lm$coefficients[1] # intercept

      # Drift corrected Mass40
      lm <- lm(calslope$calslope_40[groupNo:(groupNo+1)] ~
                 c(group.block$Time[1], tail(group.block$Time, n = 1)))
      calslope$DRIFT.calslope_40[groupNo] <- lm$coefficients[2] # slope
      calslope$DRIFT.calintercept_40[groupNo] <- lm$coefficients[1] # intercept

    }

    # 7. Interpolate calibration slopes ----------------------------

    # add columns to data
    data$INTERPOLATED.calslope_28 <- NA
    data$INTERPOLATED.calintercept_28 <- NA
    data$INTERPOLATED.calslope_32 <- NA
    data$INTERPOLATED.calintercept_32 <- NA
    data$INTERPOLATED.calslope_40 <- NA
    data$INTERPOLATED.calintercept_40 <- NA

    # create list to fill with interpolated values
    datalist = list()

    # interpolate values based on
    # calibration slope + (drift calslope *(start sample group time - sample time))
    for (groupNo in 1:max(data$Group)){

      # cut data into groups of samples
      sample.block <-
        data %>%
        filter(Type == "Sample" && Group == groupNo)

      for (i in 1:nrow(sample.block)){
        # Mass28
        sample.block$INTERPOLATED.calslope_28[i] <-
          calslope$calslope_28[groupNo] +
          (calslope$DRIFT.calslope_28[groupNo] *
              as.numeric(difftime(sample.block$Time[i],
                                  sample.block$Time[1], units = "secs")))

        sample.block$INTERPOLATED.calintercept_28[i] <-
          calslope$calintercept_28[groupNo] +
          (calslope$DRIFT.calintercept_28[groupNo] *
             as.numeric(difftime(sample.block$Time[i],
                                 sample.block$Time[1], units = "secs")))
        # Mass32
        sample.block$INTERPOLATED.calslope_32[i] <-
          calslope$calslope_32[groupNo] +
          (calslope$DRIFT.calslope_32[groupNo] *
             as.numeric(difftime(sample.block$Time[i],
                                 sample.block$Time[1], units = "secs")))

        sample.block$INTERPOLATED.calintercept_32[i] <-
          calslope$calintercept_32[groupNo] +
          (calslope$DRIFT.calintercept_32[groupNo] *
             as.numeric(difftime(sample.block$Time[i],
                                 sample.block$Time[1], units = "secs")))

        # Mass40
        sample.block$INTERPOLATED.calslope_40[i] <-
          calslope$calslope_40[groupNo] +
          (calslope$DRIFT.calslope_40[groupNo] *
             as.numeric(difftime(sample.block$Time[i],
                                 sample.block$Time[1], units = "secs")))

        sample.block$INTERPOLATED.calintercept_40[i] <-
          calslope$calintercept_40[groupNo] +
          (calslope$DRIFT.calintercept_40[groupNo] *
             as.numeric(difftime(sample.block$Time[i],
                                 sample.block$Time[1], units = "secs")))
      }

      # create list of sample blocks
      datalist[[groupNo]] <- sample.block
    }

    # convert datalist from list to dataframe
    # this dataframe will become the "detailed" data output to the user
    sampledata <- dplyr::bind_rows(datalist)

    # 8. Calculate drift and temperature corrected calibration factors -------

    # add columns for interpolated calfactors
    sampledata$INTERPOLATED.calfactor_28 <- NA
    sampledata$INTERPOLATED.calfactor_32 <- NA
    sampledata$INTERPOLATED.calfactor_40 <- NA

    # (interpolated calslope * temperature at collection) + interpolated calintercept
    # Mass 28
    sampledata$INTERPOLATED.calfactor_28 <-
      (sampledata$INTERPOLATED.calslope_28 * sampledata$CollectionTemp) +
      sampledata$INTERPOLATED.calintercept_28

    # Mass 32
    sampledata$INTERPOLATED.calfactor_32 <-
      (sampledata$INTERPOLATED.calslope_32 * sampledata$CollectionTemp) +
      sampledata$INTERPOLATED.calintercept_32

    # Mass 40
    sampledata$INTERPOLATED.calfactor_40 <-
      (sampledata$INTERPOLATED.calslope_40 * sampledata$CollectionTemp) +
      sampledata$INTERPOLATED.calintercept_40

    # 9. Calculate final concentrations -------------------------------------

    # add columns for final concentrations
    sampledata$N2conc_uM <- NA
    sampledata$O2conc_uM <- NA
    sampledata$Arconc_uM <- NA

    #BG corrected reading * interpolated calfactor
    sampledata$N2conc_uM <- sampledata$X28 * sampledata$INTERPOLATED.calfactor_28
    sampledata$O2conc_uM <- sampledata$X32 * sampledata$INTERPOLATED.calfactor_32
    sampledata$Arconc_uM <- sampledata$X40 * sampledata$INTERPOLATED.calfactor_40








    #################################################
    # Calculate concentration ratios using microM data
    data$Ratio_N2.Ar <- data$N2_microM / data$Ar_microM
    data$Ratio_O2.Ar <- data$O2_microM / data$Ar_microM

    # Calculate concentrations in mg
    data$N2_mg <- data$N2_microM * 10^(-6) * 28 * 10^3
    data$O2_mg <- data$O2_microM * 10^(-6) * 32 * 10^3
    data$Ar_mg <- data$Ar_microM * 10^(-6) * 40 * 10^3

    # Results dataframe
    # Create a dataframe of "results" that doesn't include the raw data
    calculations <-
      data[, -which(names(data) %in%
                      c("Time", "Index","X28", "X32",
                        "X40", "X99", "N2.Ar", "O2.Ar",
                        "CalFactor_28", "CalFactor_32", "CalFactor_40"))]

    # UPDATEFLAG: check the units on everything, make sure you include units in output
    output <- list(calculations = calculations,
                   # parameters = dataframe of solubility conc and cal factors
                   solubilityConcentrations = solubility.conc,
                   calibrationFactors_uM.signal = calfactor,
                   calibrationSlope = calslope,
                   all.data = data)

    return(output)

  }

  # Single point calibration --------------------------------------
  if (all(data[1:6,1] == c("Standard","Standard","Standard",
                           "Sample","Sample","Sample"))){
    stop("Single-point temperature calibration not yet supported")
    }
}
