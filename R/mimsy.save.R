#' \code{mimsy.save} Save output to an Excel workbook
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
#' mimsy.save(data, file = "results.xlsx")
#' }
#'
#' @importFrom openxlsx "write.xlsx"
#'
#' @export

mimsy.save <- function(x, file){
  # create workbook
  wb <- openxlsx::createWorkbook()

  # add worksheets
  openxlsx::addWorksheet(wb, sheetName = "Results summary")
  openxlsx::addWorksheet(wb, sheetName = "Full results")
  openxlsx::addWorksheet(wb, sheetName = "Solubility concentrations")
  openxlsx::addWorksheet(wb, sheetName = "Calibration factors")
  openxlsx::addWorksheet(wb, sheetName = "Drift corrected factors")

  # write data
  openxlsx::writeData(wb, sheet = "Results summary", x = x$results)
  openxlsx::writeData(wb, sheet = "Full results", x = x$results.full)
  openxlsx::writeData(wb, sheet = "Solubility concentrations", x = x$solubility.Concentrations)
  openxlsx::writeData(wb, sheet = "Calibration factors", x = x$calibration.Factors)
  openxlsx::writeData(wb, sheet = "Drift corrected factors", x = x$calibration.DriftCorrection)

  # save workbook
  openxlsx::saveWorkbook(wb, file = file)
}
