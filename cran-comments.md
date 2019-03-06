## Resubmission
This is a resubmission. In this version I have:

* Example file ("data.csv") now included in inst/extdata/
* Examples in mimsy.R and mimsy.save.R are no longer wrapped in \dontrun{} and are now testable
* Examples and vignettes now write files to the temporary directory using tempdir()

## Test environments
* local Windows 10 install, R 3.5.2
* ubuntu 14.04.5 (on travis-ci), R 3.5.2

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.