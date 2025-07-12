# [ATLAS9](@id YBCATLAS9)

This submodule enables interaction with the "YBC" bolometric correction (BC) grid calculated from the ATLAS9 model atmospheres of [Castelli2003](@citet), specifically the "ODFNEW" atmospheres. Details on the BC calculations are given in section 3.1 of [Chen2019](@citet). A webpage hosting and describing the atmospheres is maintained by StSCI [here](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas) and SVO [here](https://svo2.cab.inta-csic.es/theory/newov2/index.php). Robert Kurucz's webpage hosting the grid is [here](http://kurucz.harvard.edu/grids.html). These atmosphere models have been widely used since their release as they cover a wide range of stellar parameters and provide good fits to observed spectra of hotter stars (see, e.g., [Bertone2004](@citet)). Other atmospheres like [PHOENIX](@ref YBCPHOENIX) and MARCS may be preferred for lower temperature stars (for example, [Chen2019](@citet) use PHOENIX for ``T_e < 5500`` K and ATLAS9 for ``T_e > 6500`` K, interpolating between them in the overlapping region). The full range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 3,467 K | 50,118 K |
| logg   | 0.0  | 5.0   |
| \[Fe/H\] | -2.5 dex   | 0.5 dex  |
| Av     | 0.0 mag    | 20.0 mag   |
| Rv     | 3.1    | 3.1   | 

The full grid of unique values for the dependent variables is available in `BolometricCorrections.YBC.ATLAS9.gridinfo`.

```@example ybcatlas9
using BolometricCorrections
keys(BolometricCorrections.YBC.ATLAS9.gridinfo)
```

## Types

```@docs
ATLAS9YBCGrid
ATLAS9YBCTable
```

See [here](@ref ybc_systems) for a list of the photometric systems available that can be used as input arguments to [`ATLAS9YBCGrid`](@ref) and [`ATLAS9YBCTable`](@ref).

## Chemistry
These models use solar chemical abundances from [Grevesse1998](@citet). These results do not consider metal diffusion, so the photospheric solar abundances are the same as the primordial solar abundances (i.e., `X(ATLAS9Chemistry) == X_phot(ATLAS9Chemistry)`). All metallicities have scaled-solar abundance patterns.

```@docs
BolometricCorrections.YBC.ATLAS9.ATLAS9Chemistry
```

```@example ybcatlas9
chem = chemistry(ATLAS9YBCTable)
(X = X(chem), Y = Y(chem), Z = Z(chem))
```

## ATLAS9 References
This page cites the following references:

```@bibliography
Pages = ["atlas9.md"]
Canonical = false
```