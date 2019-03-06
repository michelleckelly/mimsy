## Resubmission
This is a resubmission. In this version I have:

* Example file ("data.csv") now included in inst/extdata/
* Examples in mimsy.R and mimsy.save.R are no longer wrapped in \dontrun{} and are now testable
* Examples and vignettes now write files to the temporary directory using tempdir()

############### comments from submission

Please ensure that your functions do not write by default or in your
examples/vignettes/tests in the user's home filespace. That is not allow
by CRAN policies. Please only write/save files if the user has specified
a directory. In your examples/vignettes/tests you can write to
tempdir(). E.g.

mimsy.save(data, file = file.path(tempdir(), "results.xlsx"))

Please fix and resubmit.
###############

## Test environments
* local Windows 10 install, R 3.4.1
* ubuntu 14.04.5 (on travis-ci), R 3.5.2

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.