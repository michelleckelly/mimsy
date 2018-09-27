#' \code{mimsy} Calculate dissolved gas concentrations
#'
#' Return dissolved gas concentrations in units of micromole and milligram from membrane inlet mass spectrometer (MIMS) signal data
#'
#' @param file the name of the file which the data are to be read from. File must have 'comma seperated value' ('.csv) format and follow the formatting guidelines on the Github page. If the file does not contain an absolute path, the file name is relative to the working directory.
#' @param baromet.press the ambient barometric pressure while samples processed on the MIMS. Can be a vector, if more than one reading was taken.
#' @param units the units of barometric pressure. Must be one of "atm", "hPa", "psi", "bar", or "Torr".
#' @param std.temps a numeric vector (maximum length 2) of the water temperatures (degC) of the standard baths.
#' @param bg.correct If `FALSE` (default), no background correction is applied. If `TRUE`, background correction is applied.
#' @param salinity the salinity of standards, in units of per mille. Defaults to 0.
#' @param tz a character string that specifies which time zone to parse the date with. Defaults to the user's current time zone setting. The string must be a time zone that is recognized by the user's OS.
#'
#' @return a list, containing dataframes of calculated gas concentrations, solubility concentrations, and calibration factors.
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
#' \dontrun{
#' dat <- mimsy(file = 'MIMS_data.csv', baromet.press = 981.2, units = 'hPa', std.temps = c(12.5, 15.2))
#' }
#'
#' @importFrom lubridate "mdy_hms"
#' @importFrom dplyr "group_by"
#' @importFrom dplyr "bind_rows"
#' @importFrom magrittr "%>%"
#' @importFrom stats "filter"
#' @importFrom utils "read.csv"
#'
#' @export

mimsy <- function(file, bg.correct = FALSE, baromet.press, tz = Sys.timezone(), units, salinity = 0,
    std.temps) {

    # 1. Raw data import -------------------------------------------------------------------------------
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)

    # checks check csv has correct column names check baromet.press has acceptable units check for consistent
    # standard temperatures

    # Format time column -------------------------------------------------------------------------------

    data$Time <- lubridate::mdy_hms(paste(data$RunDate, data$Time), tz = Sys.timezone())
    data$RunDate <- NULL

    # 2. Background corrections ------------------------------------------------------------------------

    if (bg.correct == FALSE) {
        # set background correction value to zero
        bg.correct <- 0
    }
    if (bg.correct != FALSE) {
        # UPDATEFLAG
        message("Background correction not yet supported.
            In meantime, please set bg.correct to FALSE")
    }

    # Barometric pressure conversion -------------------------------------------------------------------

    if (units == "atm") {
        baromet.press.atm <- mean(baromet.press)
    }
    # hPa to atm
    if (units == "hPa") {
        baromet.press.atm <- mean(baromet.press) * 0.00098692316931427
    }
    # Torr to atm
    if (units == "Torr") {
        baromet.press.atm <- mean(baromet.press) * 760
    }
    # psi to atm
    if (units == "psi") {
        baromet.press.atm <- mean(baromet.press) * 14.6959487755142
    }
    # bar to atm
    if (units == "bar") {
        baromet.press.atm <- mean(baromet.press) * 1.01325
    }
    # stop message for non-sanctioned units
    if (!(units %in% c("atm", "hPa", "Torr", "psi", "bar"))) {
        stop("Please report barometric pressure in units of `atm`, `hPa`, `psi`, `bar`, or `Torr`.")
    }

    # 3. Calculate solubilites of dissolved gas --------------------------------------------------------

    # initialize vector to store concentration values
    solubility.conc <- data.frame(O2.conc_uMol.kg = numeric(length = 2),
                                  N2.conc_uMol.kg = numeric(length = 2),
                                  Ar.conc_uMol.kg = numeric(length = 2),
                                  row.names = std.temps)

    # O2 saturation calculation ------------------------------------------------------------------------
    for (i in seq_along(std.temps)) {

        t <- std.temps[i]

        # Vapor pressure correction use the Antoine equation to calculate vapor pressure of water [bar] See
        # NIST Chemistry WebBook for general tables, these parameters valid for temperatures between -18 to
        # 100C (Stull 1947)
        vapor.press <- exp(4.6543 - (1435.264/((t + 273.15) + -64.848)))
        vapor.press <- vapor.press * 0.98692  # conversion from [bar] to [atm]

        # pressure correction [atm] = (current pressure - vapor pressure) / (standard pressure [atm] - vapor
        # pressure)
        press.corr <- (baromet.press.atm - vapor.press)/(1 - vapor.press)

        # O2 saturation calculation Combined fit coefficients [umol/kg] (Garcia and Gordon 1992, Table 1)
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
        TS <- log((298.15 - t)/(273.15 + t))  # log() == natural log (ln)
        # Salinity [per mille]
        S <- salinity

        # Calculate O2 saturation concentration at temperature and salinity (Garcia and Gordon 1992, eqn. 8)
        lnO2.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^2 + A3 * TS^3 + A4 * TS^4 + A5 * TS^5 + S * (B0 +
            B1 * TS + B2 * TS^2 + B3 * TS^3) + C0 * S^2
        O2.sat <- exp(lnO2.sat)

        # Correct O2 saturation with pressure correction, solubility.conc units [umol/kg]
        solubility.conc$O2.conc_uMol.kg[i] <- O2.sat * press.corr
    }

    # N2 saturation calculation ------------------------------------------------------------------------
    for (i in seq_along(std.temps)) {

        t <- std.temps[i]

        # Vapor pressure correction use the Antoine equation to calculate vapor pressure of water [bar] See
        # NIST Chemistry WebBook for general tables, these parameters valid for temperatures between -18 to
        # 100C (Stull 1947)
        vapor.press <- exp(4.6543 - (1435.264/((t + 273.15) + -64.848)))
        vapor.press <- vapor.press * 0.98692  # conversion from [bar] to [atm]

        # pressure correction [atm] = (current pressure - vapor pressure) /
        # (standard pressure [atm] - vapor pressure)
        press.corr <- (baromet.press.atm - vapor.press)/(1 - vapor.press)

        # N2 saturation calculation Coefficients [umol/kg]
        # (Hamme and Emerson 2004, Table 4)
        A0 <- 6.42931
        A1 <- 2.92704
        A2 <- 4.32531
        A3 <- 4.69149
        B0 <- -7.44129 * 10^-3
        B1 <- -8.02566 * 10^-3
        B2 <- -1.46775 * 10^-2
        # Scaled temperature (Hamme and Emerson 2004, eqn. 2,
        # but identical to Garcia and Gordon 1992, eqn. 8)
        TS <- log((298.15 - t)/(273.15 + t))

        # Salinity [per mille]
        S <- salinity

        # Calculate saturation concentration at temperature and salinity
        #(Hamme and erson 2004, eqn. 1)
        lnN2.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^3 + S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3)
        N2.sat <- exp(lnN2.sat)

        # Correct saturation with pressure correction, solubility.conc units are [umol/kg]
        solubility.conc$N2.conc_uMol.kg[i] <- N2.sat * press.corr
    }

    # Ar saturation calculation ----------------------------------------------------------------------
    for (i in seq_along(std.temps)) {

        t <- std.temps[i]

        # Vapor pressure correction use the Antoine equation to calculate vapor pressure of water [bar] See
        # NIST Chemistry WebBook for general tables, these parameters valid for temperatures between -18 to
        # 100C (Stull 1947)
        vapor.press <- exp(4.6543 - (1435.264/((t + 273.15) + -64.848)))
        vapor.press <- vapor.press * 0.98692  # conversion from [bar] to [atm]

        # pressure correction [atm] = (current pressure - vapor pressure) /
        # (standard pressure [atm] - vapor pressure)
        press.corr <- (baromet.press.atm - vapor.press)/(1 - vapor.press)

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
        TS <- log((298.15 - t)/(273.15 + t))

        # Salinity [per mille]
        S <- salinity

        # Calculate saturation concentration at temperature and salinity
        # (Hamme and erson 2004, eqn. 1)
        lnAr.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^3 + S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3)
        Ar.sat <- exp(lnAr.sat)
        # Correct saturation with pressure correction, solubility.conc units are [umol/kg]
        solubility.conc$Ar.conc_uMol.kg[i] <- Ar.sat * press.corr
    }

    # 4. Calculate calibration factors ----------------------------------------------------------------

    # Group data by Type (`Standard` or `Sample`) and Group (numeric, 1:n)
    data <- dplyr::group_by(data, data$Type, data$Group)

    # two-point calibration has 6 standard readings
    if (all(data[1:6, 1] == c("Standard", "Standard", "Standard", "Standard", "Standard", "Standard"))) {

        message("Calculating dissolved concentrations based on a two-point temperature standard.")

        # assemble empty data frame for calibration factors
        calfactor <- data.frame(calfactor_28 = numeric(length = max(data$Group) * 2),
                                calfactor_32 = numeric(length = max(data$Group) *2),
                                calfactor_40 = numeric(length = max(data$Group) * 2),
                                row.names = paste0(std.temps, "degC_",
                                                   "Group_", rep(1:max(data$Group), each = 2)))

        for (groupNo in 1:max(data$Group)) {
            # individually extract each group of standards
            cal.block <- data %>% filter(data$Type == "Standard" && data$Group == groupNo)

            # calculate calibration factor calfactor = solubility concentration at std temp / avg(MIMS readings at
            # std temp)

            # Mass28 (N2) Standard temp 1
            calfactor$calfactor_28[2 * groupNo - 1] <-
              solubility.conc$N2.conc_uMol.kg[1]/mean(cal.block$X28[1:3])
            # standard temp 2
            calfactor$calfactor_28[2 * groupNo] <-
              solubility.conc$N2.conc_uMol.kg[2]/mean(cal.block$X28[4:6])

            # Mass32 (O2)
            calfactor$calfactor_32[2 * groupNo - 1] <-
              solubility.conc$O2.conc_uMol.kg[1]/mean(cal.block$X32[1:3])
            calfactor$calfactor_32[2 * groupNo] <-
              solubility.conc$O2.conc_uMol.kg[2]/mean(cal.block$X32[4:6])

            # Mass40 (Ar)
            calfactor$calfactor_40[2 * groupNo - 1] <-
              solubility.conc$Ar.conc_uMol.kg[1]/mean(cal.block$X40[1:3])
            calfactor$calfactor_40[2 * groupNo] <-
              solubility.conc$Ar.conc_uMol.kg[2]/mean(cal.block$X40[4:6])
        }

        # 5. Calculate slope and intercepts of calibration curve -----------------
        calslope <- data.frame(calslope_28 = numeric(length = max(data$Group)),
                               calintercept_28 = numeric(length = max(data$Group)),
                               calslope_32 = numeric(length = max(data$Group)),
                               calintercept_32 = numeric(length = max(data$Group)),
                               calslope_40 = numeric(length = max(data$Group)),
                               calintercept_40 = numeric(length = max(data$Group)),
                               row.names = paste0("Group", 1:max(data$Group)))

        for (groupNo in 1:max(data$Group)) {
            # using a linear model to find slope and intercept of calibration line linear model: y ~ x

            # Mass28
            lm <- lm(c(calfactor$calfactor_28[2 * groupNo - 1], calfactor$calfactor_28[2 * groupNo]) ~
                std.temps)
            calslope$calslope_28[groupNo] <- lm$coefficients[2]  # slope
            calslope$calintercept_28[groupNo] <- lm$coefficients[1]  #intercept

            # Mass32
            lm <- lm(c(calfactor$calfactor_32[2 * groupNo - 1], calfactor$calfactor_32[2 * groupNo]) ~
                std.temps)
            calslope$calslope_32[groupNo] <- lm$coefficients[2]  # slope
            calslope$calintercept_32[groupNo] <- lm$coefficients[1]  #intercept

            # Mass40
            lm <- lm(c(calfactor$calfactor_40[2 * groupNo - 1], calfactor$calfactor_40[2 * groupNo]) ~
                std.temps)
            calslope$calslope_40[groupNo] <- lm$coefficients[2]  # slope
            calslope$calintercept_40[groupNo] <- lm$coefficients[1]  #intercept
        }

        # 6. Perform drift correction for calibration slope and intercept --------

        # add columns to calslope
        calslope$DRIFT.calslope_28 <- NA
        calslope$DRIFT.calintercept_28 <- NA
        calslope$DRIFT.calslope_32 <- NA
        calslope$DRIFT.calintercept_32 <- NA
        calslope$DRIFT.calslope_40 <- NA
        calslope$DRIFT.calintercept_40 <- NA

        for (groupNo in 1:max(data$Group)) {

            # take the slope between successive calibration (slope or intercept) values

            # Drift corrected Mass28 slope
            calslope$DRIFT.calslope_28[groupNo] <- (calslope$calslope_28[groupNo + 1] - calslope$calslope_28[groupNo])/as.numeric(difftime(data$Time[data$Group ==
                (groupNo + 1)][1], data$Time[data$Group == groupNo][1], units = "days"))
            # intercept
            calslope$DRIFT.calintercept_28[groupNo] <- (calslope$calintercept_28[groupNo + 1] - calslope$calintercept_28[groupNo])/as.numeric(difftime(data$Time[data$Group ==
                (groupNo + 1)][1], data$Time[data$Group == groupNo][1], units = "days"))

            # Drift corrected Mass32 slope
            calslope$DRIFT.calslope_32[groupNo] <- (calslope$calslope_32[groupNo + 1] - calslope$calslope_32[groupNo])/as.numeric(difftime(data$Time[data$Group ==
                (groupNo + 1)][1], data$Time[data$Group == groupNo][1], units = "days"))
            # intercept
            calslope$DRIFT.calintercept_32[groupNo] <- (calslope$calintercept_32[groupNo + 1] - calslope$calintercept_32[groupNo])/as.numeric(difftime(data$Time[data$Group ==
                (groupNo + 1)][1], data$Time[data$Group == groupNo][1], units = "days"))

            # Drift corrected Mass40 slope
            calslope$DRIFT.calslope_40[groupNo] <- (calslope$calslope_40[groupNo + 1] - calslope$calslope_40[groupNo])/as.numeric(difftime(data$Time[data$Group ==
                (groupNo + 1)][1], data$Time[data$Group == groupNo][1], units = "days"))
            # intercept
            calslope$DRIFT.calintercept_40[groupNo] <- (calslope$calintercept_40[groupNo + 1] - calslope$calintercept_40[groupNo])/as.numeric(difftime(data$Time[data$Group ==
                (groupNo + 1)][1], data$Time[data$Group == groupNo][1], units = "days"))

            # if there's no standard group at the tail (aka, run ends after a sample) base the drift correction on
            # the preceding block of standards
            if (!is.element(groupNo + 1, unique(data$Group))) {
                # copy down drift correction from preceeding group this isn't the best fix in the world, but better
                # than NA for now
                calslope$DRIFT.calslope_28[groupNo] <- calslope$DRIFT.calslope_28[groupNo - 1]
                calslope$DRIFT.calintercept_28[groupNo] <- calslope$DRIFT.calintercept_28[groupNo - 1]
                calslope$DRIFT.calslope_32[groupNo] <- calslope$DRIFT.calslope_32[groupNo - 1]
                calslope$DRIFT.calintercept_32[groupNo] <- calslope$DRIFT.calintercept_32[groupNo - 1]
                calslope$DRIFT.calslope_40[groupNo] <- calslope$DRIFT.calslope_40[groupNo - 1]
                calslope$DRIFT.calintercept_40[groupNo] <- calslope$DRIFT.calintercept_40[groupNo - 1]
            }
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

        # interpolate values based on calibration slope + (drift calslope *(start sample group time - sample
        # time))
        for (groupNo in 1:max(data$Group)) {

            # cut data into groups of samples
            group.block <- data %>% filter(data$Group == groupNo)

            for (i in 1:nrow(group.block)) {
                # Mass28
                group.block$INTERPOLATED.calslope_28[i] <- calslope$calslope_28[groupNo] + (calslope$DRIFT.calslope_28[groupNo] *
                  as.numeric(difftime(group.block$Time[i], group.block$Time[1], units = "days")))

                group.block$INTERPOLATED.calintercept_28[i] <- calslope$calintercept_28[groupNo] + (calslope$DRIFT.calintercept_28[groupNo] *
                  as.numeric(difftime(group.block$Time[i], group.block$Time[1], units = "days")))
                # Mass32
                group.block$INTERPOLATED.calslope_32[i] <- calslope$calslope_32[groupNo] + (calslope$DRIFT.calslope_32[groupNo] *
                  as.numeric(difftime(group.block$Time[i], group.block$Time[1], units = "days")))

                group.block$INTERPOLATED.calintercept_32[i] <- calslope$calintercept_32[groupNo] + (calslope$DRIFT.calintercept_32[groupNo] *
                  as.numeric(difftime(group.block$Time[i], group.block$Time[1], units = "days")))

                # Mass40
                group.block$INTERPOLATED.calslope_40[i] <- calslope$calslope_40[groupNo] + (calslope$DRIFT.calslope_40[groupNo] *
                  as.numeric(difftime(group.block$Time[i], group.block$Time[1], units = "days")))

                group.block$INTERPOLATED.calintercept_40[i] <- calslope$calintercept_40[groupNo] + (calslope$DRIFT.calintercept_40[groupNo] *
                  as.numeric(difftime(group.block$Time[i], group.block$Time[1], units = "days")))
            }  # close internal row for loop

            # create list of sample blocks
            datalist[[groupNo]] <- group.block

        }  # close external group for loop

        # convert datalist from list to dataframe this dataframe will become the 'detailed' data output to the
        # user
        data <- dplyr::bind_rows(datalist)

        # 8. Calculate drift and temperature corrected calibration factors -------

        # (interpolated calslope * temperature at collection) + interpolated calintercept Mass 28
        data$INTERPOLATED.calfactor_28 <- (data$INTERPOLATED.calslope_28 * data$CollectionTemp) + data$INTERPOLATED.calintercept_28

        # Mass 32
        data$INTERPOLATED.calfactor_32 <- (data$INTERPOLATED.calslope_32 * data$CollectionTemp) + data$INTERPOLATED.calintercept_32

        # Mass 40
        data$INTERPOLATED.calfactor_40 <- (data$INTERPOLATED.calslope_40 * data$CollectionTemp) + data$INTERPOLATED.calintercept_40

        # 9. Calculate final concentrations -------------------------------------

        # BG corrected reading * interpolated calfactor
        data$N2_uMol <- data$X28 * data$INTERPOLATED.calfactor_28
        data$O2_uMol <- data$X32 * data$INTERPOLATED.calfactor_32
        data$Ar_uMol <- data$X40 * data$INTERPOLATED.calfactor_40

        # convert from microM to mg
        data$N2_mg <- data$N2_uMol * 10^(-6) * 28 * 10^3
        data$O2_mg <- data$O2_uMol * 10^(-6) * 32 * 10^3
        data$Ar_mg <- data$Ar_uMol * 10^(-6) * 40 * 10^3

        # 10. Output results to user -------------------------------------------

        # results A subset of `data` that contains only sample identifiers and final concentration results.
        # This is a more user-friendly short-form

        # i'm subtracting out names rather than selecting for names, because I want to preserve any other
        # sample identity columns that the user may have added to the orignal .csv
        results <- data[, -which(names(data) %in% c("Index", "Time", "X28", "X32", "X40", "X99", "N2.Ar",
            "O2.Ar", "INTERPOLATED.calslope_28", "INTERPOLATED.calintercept_28", "INTERPOLATED.calslope_32",
            "INTERPOLATED.calintercept_32", "INTERPOLATED.calslope_40", "INTERPOLATED.calintercept_40",
            "INTERPOLATED.calfactor_28", "INTERPOLATED.calfactor_32", "INTERPOLATED.calfactor_40"))]
        # grab only the sample results
        results <- results[results$Type == "Sample", ]
        results$Type <- NULL
        results$Group <- NULL

        outlist <- list(results = results, solubility.Concentrations = solubility.conc, calibration.Factors = calfactor,
            calibration.DriftCorrection = calslope, results.full = data)
        return(outlist)
    }
    # if (all(data[1:6,1] == c('Standard','Standard','Standard','Sample','Sample','Sample'))) {
    # stop('Single-point temperature calibration not yet supported') }
}
