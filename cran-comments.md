## Test environments
* local Windows 10 Pro Install, R 3.4.1

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:
* checking R code for possible problems ... NOTE
  mimsy: no visible binding for global variable 'Type'
  mimsy: no visible binding for global variable 'Group'
  Undefined global functions or variables:
    Group Type
	
Global variables "Type" and "Group" are defined in the specially formatted CSV 
input file that the user must upload. A detailed description of this formatting
is included in /vignettes/mimsy.Rmd