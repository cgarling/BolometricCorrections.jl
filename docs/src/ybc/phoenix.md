# [PHOENIX](@id YBCPHOENIX)

This submodule enables interaction with the "YBC" bolometric correction (BC) grid of PHOENIX model atmospheres described in section 3.2 of [Chen2019](@citet). The original atmosphere models are hosted by SVO [here](https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=bt-settl) and StSCI provides a subset of the PHOENIX library for use with their Synphot software [here](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot).

In the original YBC implementation, the PHOENIX models are used for cool stars ``T_e < 5500`` K and transition smoothly onto the [ATLAS9](@ref YBCATLAS9) models from ``5500 < T_e < 6500`` K. YBC specifically uses the BT-Settl set of PHOENIX models for which the main reference is [Allard2012](@citet). The coverage of the model grid is fairly large -- the full range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 2570 K | 70,794 K |
| logg   | -0.5   | 6.0   |
| \[Fe/H\] | -4.0 dex   | 0.5 dex  |
| Av     | 0.0 mag    | 20.0 mag   |
| Rv     | 3.1    | 3.1   | 

The full grid of unique values for the dependent variables is available in `BolometricCorrections.MIST.gridinfo`.

```@example ybcphoenix
using BolometricCorrections
keys(BolometricCorrections.YBC.PHOENIX.gridinfo)
```

## Types

```@docs
PHOENIXYBCGrid
PHOENIXYBCTable
```

See [here](@ref ybc_systems) for a list of the photometric systems available that can be used as input arguments to [`PHOENIXYBCGrid`](@ref) and [`PHOENIXYBCTable`](@ref). 

## Chemistry
These models use solar chemical abundances from [Asplund2009](@citet), the same as [MIST](@ref MIST), so we use the same chemistry type to represent the PHOENIX chemistry (see [here](@ref MIST_chemistry) for more information).

```@example ybcphoenix
chemistry(PHOENIXYBCTable)
```

[Allard2013](@citet) discusses alternative BT-Settl models with the solar chemical composition from [Caffau2011](@citet), YBC uses the models with [Asplund2009](@citet) abundances.


## PHOENIX References
This page cites the following references:

```@bibliography
Pages = ["phoenix.md"]
Canonical = false
```