#' \code{mimsy.save} Save MIMS calculations to an Excel workbook
#'
#' Easily save the full output of the `mimsy` calculation function to a nicely-formatted, multi-tab Excel .xlsx file.
#'
#' @param x object to be written to file, the output of `mimsy` function
#' @param file desired file name with .xlsx ending. Example: "mimsyCalculations.xlsx"
#'
#' @return Outputs an excel workbook to the user's working directory
#'
#' @examples
#' \dontrun{
#' mimsy.save(dat, file = "mimsyCalculations.xlsx")
#' }
#'
#' @importFrom openxlsx "write.xlsx"
#'
#' @export
mimsy.save <- function(x, file){
  # seperate list output into individual dataframes
  results <- as.data.frame(x$results)
  fullresults <- as.data.frame(x$results.full)
  solcon <- as.data.frame(x$solubility.Concentrations)
  calfac <- as.data.frame(x$calibration.Factors)
  driftcalfac <- as.data.frame(x$calibration.DriftCorrection)

  openxlsx::write.xlsx(results, file = file, sheetName = "Results summary",
                       col.names = TRUE, row.names = FALSE, append = FALSE)
  openxlsx::write.xlsx(fullresults, file = file,
                       sheetName = "Full results", col.names = TRUE,
                       row.names = FALSE, append = TRUE)
  openxlsx::write.xlsx(solcon, file = file,
                       sheetName = "Solubility concentrations", col.names = TRUE,
                       row.names = TRUE, append = TRUE)
  openxlsx::write.xlsx(calfac, file = file,
                       sheetName = "Calibration factors", col.names = TRUE,
                       row.names = TRUE, append = TRUE)
  openxlsx::write.xlsx(driftcalfac, file = file,
                       sheetName = "Drift corrected calibr factors", col.names = TRUE,
                       row.names = TRUE, append = TRUE)
}
