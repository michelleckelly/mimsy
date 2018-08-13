# mimsy
## The unofficial MIMS data cruncher

### Description
`mimsy` is the unofficial R package for crunching raw MIMS (Membrane Inlet Mass Spectrometer, 
Bay Instruments) day. Use the `mimsy` function to calculate dissolved gas concentrations and 
ratios (N2, O2, Ar) based on standard gas solubility at specific temperature, pressure, and salinity. 
Supports dual-temperature standard calibration, for setups with two water baths, and outputs dissolved 
gas concentration in units of microM or mg.

### Installation
```R
library(devtools)
install_github("michelleckelly/mimsy")
```

### Example

### Working task list
- [x] saturation calculations for O2, N2, Ar scripted and references cited
- [ ] settle on indexing choice and script calcibration factors and slope corrections
- [ ] add additional unit conversion options for barometric pressure
- [ ] add support for background corrections
- [ ] add support for CH4
- [ ] write input .csv formatting guidelines
- [ ] write short example for README
- [ ] write full vignette
