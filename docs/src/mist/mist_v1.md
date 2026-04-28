# [MIST v1.2](@id MIST_v1)

The MIST v1.2 BC grid covers a wide range of photometric filters and spans the full range of effective temperature, surface gravity, and metallicity relevant for most stellar evolution applications. The following figure shows a projection of a small portion of the BC table for one choice of metallicity and V-band extinction.

```@example
using BolometricCorrections # hide
include(joinpath(@__DIR__, "..", "plots.jl")) # hide
grid = MISTv1BCGrid("JWST") # hide
Teff = logrange(exp10(3.5), 10_000; length=1000) # hide
# logg = range(extrema(MISTv1BCGrid).logg...; length=1000) # hide
logg = range(3.0, 5.0; length=1000) # hide
f, ax = plot_bc_table(grid(-1, 0), "F090W", Teff, logg) # hide
ax.title = "MIST v1.2 BCs for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text="[M/H] = -1\n Av = 0", align=(:right, :top), space=:relative) # hide
f # hide
```

## [Chemistry](@id MIST_chemistry)

The MIST v1.2 BC grid assumes scaled-solar metal abundance ratios using the protostellar solar abundances of [Asplund2009](@citet), so \[M/H\] is equivalent to \[Fe/H\]. We provide [`BolometricCorrections.MIST.MISTv1Chemistry`](@ref) to access information on the v1.2 chemical mixture following the [chemical mixture API](@ref chemistry_api). `MISTChemistry` is provided as a deprecated alias for people transitioning from earlier versions of BolometricCorrections.jl prior to the addition of MIST v2.5.

```@docs
BolometricCorrections.MIST.MISTv1Chemistry
```

The full range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 2,500 K | 1,000,000 K |
| logg   | −4.0   | 9.5   |
| \[Fe/H\] | −4.0 dex   | +0.75 dex  |
| Av     | 0.0 mag    | 6.0 mag   |
| Rv     | 3.1    | 3.1   |

The full grid of unique values for the dependent variables is available in `BolometricCorrections.MIST.gridinfov1`.

```@example
import BolometricCorrections # hide
keys(BolometricCorrections.MIST.gridinfov1)
```

## Types

```@docs
BolometricCorrections.MIST.MISTv1BCGrid
```

The constructor for [`MISTv1BCGrid`](@ref) includes a parser to translate most human-readable photometric system names (e.g., `"HST/ACS-WFC"`) into their proper internal identifiers. The full list of internal specifiers is given below.

```@example
import DataDeps # hide
using BolometricCorrections.MIST # hide
# custom `show` call to prevent truncation of output # hide
show(stdout, "text/plain", filter(x -> occursin("MIST", x) & !occursin("2.5", x), keys(DataDeps.registry)))
```

Each photometric system uses a separate data file. The first time you request a BC grid for a system you have not yet downloaded, you will be prompted to allow the download. Once a grid is constructed for a particular photometric system, a BC table with fixed \[Fe/H\] and ``A_V`` can be interpolated.

```@docs
BolometricCorrections.MIST.MISTv1BCTable
```

### Deprecated aliases

`MISTBCGrid` and `MISTBCTable` are deprecated aliases for [`MISTv1BCGrid`](@ref) and [`MISTv1BCTable`](@ref), respectively, provided for backwards compatibility. They forward all arguments to the v1 constructors and emit a deprecation warning if Julia is configured to show them (e.g., Julia started as `julia --depwarn=yes`).

## Photometric Zeropoints

The MIST v1.2 bolometric corrections adopt a solar bolometric magnitude of ``M_{\odot,\text{bol}} = 4.74`` and a solar bolometric luminosity of ``L_\odot = 3.828 \times 10^{33}`` erg s⁻¹. Information needed to convert between the AB, Vega, and ST photometric systems is contained in [`BolometricCorrections.MIST.zpt`](@ref) (documented on the [MIST overview page](@ref MIST)).

## MIST v1.2 References
This page cites the following references:

```@bibliography
Pages = ["mist_v1.md"]
Canonical = false
```
