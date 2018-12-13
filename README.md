[![PyPI status](https://img.shields.io/pypi/status/ansicolortags.svg)](https://GitHub.com/michelleckelly/mimsy/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/michelleckelly/mimsy/graphs/commit-activity)

## mimsy
_calculate dissolved gas concentrations from MIMS signal data without getting a headache_

### Description  
`mimsy` is an R package that calculates dissolved gas concentrations from raw MIMS (Membrane Inlet Mass Spectrometer, Bay Instruments) signal data based on dissolved gas solubilities. Use the `mimsy` function on a formatted .csv file to return dissolved gas concentrations (mg and μMole) of N<sub>2</sub>, O<sub>2</sub>, Ar based on gas solubility at temperature, pressure, and salinity. Then, easily save the output of `mimsy` to a multi-tab Excel workbook with the `mimsy.save` function. Supports dual-temperature standard calibration for dual-bath MIMS setups.

### Installation  

```R
# Pull package from github using devtools
library(devtools)
install_github("michelleckelly/mimsy", dependencies = "Depends")

# Load package into your R environment
library(mimsy)
```

### Citation
Run `citation("mimsy")` in the R console to view the recommended citation for this package.

### Disclaimer
`mimsy` holds no official endorsement from the Bay Instruments company. This software is preliminary and subject to revision. By the use of this software, the user assumes their own responsibility for ensuring the accuracy of the program. 

### References
Garcia, H., and L. Gordon (1992), _Oxygen solubility in seawater: Better fitting
equations._ Limnology and Oceanography, 37(6). https://doi.org/10.4319/lo.1992.37.6.1307

Benson, B. B. & Krause, D. (1984). _The concentration and isotopic
fractionation of oxygen dissolved in freshwater and seawater in equilibrium
with the atmosphere._ Limnology and Oceanography, 29(3), 620-632.
https://doi.org/10.4319/lo.1984.29.3.0620

Stull, D. R. (1947). _Vapor Pressure of Pure Substances. Organic and
Inorganic Compounds._ Industrial & Engineering Chemistry, 39(4), 517-540.
https://doi.org/10.1021/ie50448a022

Hamme, R. C. & Emerson, S. R. (2004). _The solubility of neon, nitrogen and argon
in distilled water and seawater._ Deep-Sea Research I, 51(11), 1517-1528. 
https://doi.org/10.1016/j.dsr.2004.06.009
