#' \code{mimsy.save} Save full output of `mimsy` to a multi-tabbed Excel workbook
#'
#' Easily save the full output of the `mimsy` calculation function to a nicely-formatted, multi-tab Excel .xlsx file.
#'
#' @param x object to be written to file, the output of `mimsy` function
#' @param file desired file name with .xlsx ending. Example: "mimsyCalculations.xlsx"
#'
#' @return Excel workbook to the user's working directory
#'
#' @examples
#' mimsy.save(myData, file = "mimsyCalculations.xlsx")
#'
#' @importFrom xlsx "write.xlsx"
#'
#' @export

mimsy.save <- function(x, file){
  # seperate list output into individual dataframes
  results <- as.data.frame(x$results)
  fullresults <- as.data.frame(x$results.full)
  solcon <- as.data.frame(x$solubility.Concentrations)
  calfac <- as.data.frame(x$calibration.Factors)
  driftcalfac <- as.data.frame(x$calibration.DriftCorrection)

  xlsx::write.xlsx(results, file = file, sheetName = "Results summary",
                   col.names = TRUE, row.names = FALSE, append = FALSE)
  xlsx::write.xlsx(fullresults, file = file,
             sheetName = "Full results", col.names = TRUE, row.names = FALSE, append = TRUE)
  xlsx::write.xlsx(solcon, file = file,
             sheetName = "Solubility concentrations", col.names = TRUE, row.names = TRUE, append = TRUE)
  xlsx::write.xlsx(calfac, file = file,
                   sheetName = "Calibration factors", col.names = TRUE,
                   row.names = TRUE, append = TRUE)
  xlsx::write.xlsx(driftcalfac, file = file,
             sheetName = "Drift corrected calibr factors", col.names = TRUE,
             row.names = TRUE, append = TRUE)
}
