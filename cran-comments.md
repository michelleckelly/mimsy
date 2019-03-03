## Test environments
* local Windows 10 install, R 3.4.1
* ubuntu 14.04.5 (on travis-ci), R 3.5.2

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:
* checking R code for possible problems ... NOTE
  mimsy: no visible binding for global variable 'Type'
  mimsy: no visible binding for global variable 'Group'
  Undefined global functions or variables:
    Group Type
	
Global variables "Type" and "Group" are defined in the specially formatted CSV 
input file that the user must upload (as specific column names). A detailed description 
of this required formatting is included in /vignettes/mimsy.Rmd

## Downstream dependencies
There are currently no downstream dependencies for this package.