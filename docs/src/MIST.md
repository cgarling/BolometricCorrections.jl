# [MIST](@id MIST)

This submodule enables interaction with the bolometric correction (BC) grid released as part of the Mesa Isochrones & Stellar Tracks ([MIST](https://waps.cfa.harvard.edu/MIST/) [Dotter2016,Choi2016](@cite)) project. The MIST BC grid is convenient because it includes a wide range of photometric filters and covers the full range of effective temperature, surface gravity, and metallicity relevant for most applications in stellar evolution. The grid is also *regular* in the dependent variables, greatly simplifying interpolation. The following figure shows a projection of a small portion of the BC table for one choice of metallicity and V-band extinction.

```@setup mist_plotting
include(joinpath(@__DIR__, "examples", "bc_tables.jl"))
```
```@example mist_plotting
plot_mist_bc_table("JWST", "F090W", -1, 0) # hide
```

The MIST BC grid assumes scaled-solar metal abundance ratios assuming the protostellar birth cloud bulk metallicity of [Asplund2009](@citet), so \[M/H\] is equivalent to \[Fe/H\]. The literature on the MIST models prefers to use \[Fe/H\], so we follow the same convention here.

The full range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 2500 K | 1e6 K |
| logg   | -4.0   | 9.5   |
| \[Fe/H\] | -4.0 dex   | 0.75 dex  |
| Av     | 0.0 mag    | 6.0 mag   |
| Rv     | 3.1    | 3.1   |

The full grid of unique values for the dependent variables is availabie in `BolometricCorrections.MIST.gridinfo`.

```@example
import BolometricCorrections # hide
keys(BolometricCorrections.MIST.gridinfo)
```

## Types

```@docs
MISTBCGrid
```

The constructor for `MISTBCGrid` includes a parser to translate most human-readable photometric system names (e.g., `"HST/ACS-WFC"`) into their proper internal identifiers. The full list of internal specifiers is given below.

```@example
import DataDeps
using BolometricCorrections.MIST
# custom `show` call to prevent truncation of output # hide
show(stdout, "text/plain", keys(DataDeps.registry))
```

Once a BC grid has been constructed for a particular choice of photometric system, a BC table (with variables \[Fe/H\] and Av fixed) can be interpolated. 

```@docs
MISTBCTable
```

## Photometric Zeropoints
The MIST bolometric corrections assume a bolometric luminosity zeropoint of ``3.0128 \times 10^{35} \, \text{erg} \, \text{s}^{-1}`` to define ``M_\text{bol} = 0``. This is equivalent to adopting solar values for the bolometric magnitude of ``M_\text{bol} = 4.74`` mag with bolometric luminosity of ``3.828 \times 10^{33} \, \text{erg} \, \text{s}^{-1}``.

Information needed to convert between different photometric systems (AB, Vega, ST) is contained in [`BolometricCorrections.MIST.zeropoints`](@ref), which is an instance of the [`BolometricCorrections.MIST.MISTZeropoints`](@ref) type. Additional information on operations supported by this type is available in our [API documentation](@ref zpt_api).

```@docs
BolometricCorrections.MIST.zeropoints
BolometricCorrections.MIST.MISTZeropoints
```

## MIST References
This page cites the following references:

```@bibliography
Pages = ["MIST.md"]
Canonical = false
```