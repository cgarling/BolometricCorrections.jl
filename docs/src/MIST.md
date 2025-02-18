# [MIST](@id MIST)

This submodule enables interaction with the bolometric correction (BC) grid released as part of the Mesa Isochrones & Stellar Tracks ([MIST](https://waps.cfa.harvard.edu/MIST/)) project. The range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 2500 K | 1e6 K |
| logg   | -4.0   | 9.5   |
| \[Fe/H\] | -4.0   | 0.75  |
| Av     | 0.0    | 6.0   |
| Rv     | 3.1    | 3.1   |

The full grid of unique values for the dependent variables is availabie in `BolometricCorrections.MIST.gridinfo`.

```@example
import BolometricCorrections # hide
keys(BolometricCorrections.MIST.gridinfo)
```

# Types

```@docs
MISTBCGrid
MISTBCTable
```

# Photometric Zeropoints
```@docs
BolometricCorrections.MIST.zeropoints
BolometricCorrections.MIST.MISTZeropoints
```

# References
 - See section 5.4 of [Choi+2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...823..102C/abstract) for a description of the MIST bolometric corrections.