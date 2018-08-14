# `mimsy`  

### Description  

`mimsy` is the unofficial R package for transforming raw MIMS (Membrane Inlet Mass Spectrometer, Bay Instruments) data into dissolved gas concentration readings. The `mimsy` function calculate dissolved gas concentrations and ratios (N2, O2, Ar) based on gas solubility at specific temperature, pressure, and salinity. Oxygen solubility calculations based on Garcia and Gordon (1992); nitrogen and argon solubility calculations based on Hamme and Emerson (2004). See `mimsy` for full reference list. Supports dual-temperature standard calibration.

### Installation  

```R
# Pull package from github using devtools
library(devtools)
install_github("michelleckelly/mimsy")
# Load package into R environment
library(mimsy)
```

### Example  

`mimsy` is currently in development. A simple vignette will be added here once the package is ready to rumble

### Working task list  

- [x] saturation calculations for O2, N2, Ar scripted and references cited
- [ ] settle on indexing choice and script calcibration factors and slope corrections
- [ ] add additional unit conversion options for barometric pressure
- [ ] add support for background corrections
- [ ] add support for CH4
- [ ] write input .csv formatting guidelines
- [ ] write short example for README
- [ ] write full vignette
