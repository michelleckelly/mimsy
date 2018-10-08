[![PyPI status](https://img.shields.io/pypi/status/ansicolortags.svg)](https://GitHub.com/michelleckelly/mimsy/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/michelleckelly/mimsy/graphs/commit-activity)

### Description  

Calculate dissolved gas concentrations from raw MIMS (Membrane Inlet Mass Spectrometer, Bay Instruments) signal data. Use the `mimsy` function on a formatted .csv file to return dissolved gas concentrations (mg and μMole) and ratios of N<sub>2</sub>, O<sub>2</sub>, Ar based on gas solubility at temperature, pressure, and salinity. Then, easily save the output of `mimsy` to a nicely-formatted multi-tab Excel workbook with the `mimsy.save` function. Supports dual-temperature standard calibration for dual-bath MIMS setups.

### Installation  

```R
# Pull package from github using devtools
library(devtools)
install_github("michelleckelly/mimsy")
#
# Load package into your R environment
library(mimsy)
```

### References
Garcia, H., and L. Gordon (1992), _Oxygen solubility in seawater: Better fitting
equations._ Limnology and Oceanography, 37(6).

Benson, B. B. & Krause, D. (1984). _The concentration and isotopic
fractionation of oxygen dissolved in freshwater and seawater in equilibrium
with the atmosphere._ Limnology and Oceanography, 29(3), 620-632.
doi:10.4319/lo.1984.29.3.0620

Stull, D. R. (1947). _Vapor Pressure of Pure Substances. Organic and
Inorganic Compounds._ Industrial & Engineering Chemistry, 39(4), 517-540.
doi: 10.1021/ie50448a022

Hamme, R. C. & Emerson, S. R. (2004). _The solubility of neon, nitrogen and argon
in distilled water and seawater._, Deep-Sea Research I, 51(11), 1517-1528.
