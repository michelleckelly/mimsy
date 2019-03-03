
# mimsy 0.4.2

## Minor changes
* Preparation for submission to CRAN
* `README.md` updated with full workflow example
* logo updated

# mimsy 0.4.1

## Major changes
* `mimsy()` takes `data` input as a dataframe, no longer internally loads from file path

## Minor changes
* `vignettes/mimsy.Rmd` rewritten for clarity and screenshots added

# mimsy 0.3.1

## Major changes
* `mimsy()` std.temps pulled from CollectionTemps column of Group
1 Standards, std.temps input argument now depreciated

## Minor changes
* `mimsy()` now indicates the standard temperatures used in 
an output message to the user
* `mimsy()` now checks that the length of std.temps vector
is appropriate

## Bug fixes
* update `mimsy.save()` to rely on non-depreciated forms of 
`openxlsx` functions

# mimsy 0.2.1

## Minor changes
* `vignettes/mimsy.Rmd` now includes directions to save output
as RData file

## Bug fixes
* `vignettes/mimsy.Rmd` and `README.md` install instructions
updated to fix bug where depenencies weren't downloading properly.

# mimsy 0.2.0

## Major changes
* Added `vignettes/mimsy.Rmd` quick start guide for package.
* Added `NEWS.md` file to track changes to the package.
