# [MIST v2.5](@id MIST_v2)

The MIST v2.5 BC grid extends the [v1.2 grid](@ref MIST_v1) by introducing \[α/Fe\] as a free parameter and assumes a different solar chemistry than v1.2. It also expands the photometric system coverage and revises the grid of dependent variables (temperature, surface gravity, etc.). The following figure shows a projection of a small portion of the BC table for one choice of metallicity, α-enhancement, and V-band extinction.

```@example mistv2_plot
using BolometricCorrections # hide
include(joinpath(@__DIR__, "..", "plots.jl")) # hide
grid = MISTBCGridv2("JWST") # hide
feh, afe, Av = -1, 0, 0 # hide
Teff = logrange(exp10(3.5), 10_000; length=1000) # hide
# logg = range(extrema(MISTBCGridv2).logg...; length=1000) # hide
logg = range(3.0, 5.0; length=1000) # hide
table = grid(feh, afe, Av) # hide
f, ax = plot_bc_table(table, "F090W", Teff, logg) # hide
ax.title = "MIST v2.5 BCs for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text="[Fe/H] = $(feh)\n [α/Fe] = $(afe)\n Av = $(Av)", align=(:right, :top), space=:relative) # hide
f # hide
```

and here is the difference between the MIST v1.2 and v2.5 BCs for this filter and metallicity. Over this range of temperature and surface gravity, the v1.2 and v2.5 BCs are very similar.

```@example mistv2_plot
grid_v1 = MISTBCGridv1("JWST") # hide
table_v1 = grid_v1(feh, Av) # hide
f, ax = plot_bc_table_diff(table_v1, table, ("F090W", "F090W"), Teff, logg) # hide
ax.title = "MIST v1.2 - v2.5 for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text="[Fe/H] = $(feh)\n [α/Fe] = $(afe)\n Av = $(Av)", align=(:right, :top), space=:relative) # hide
f # hide
```

## [Chemistry](@id MIST_v2_chemistry)

The MIST v2.5 BC grid uses solar chemical abundances from [Grevesse1998](@citet), in contrast to [Asplund2009](@citet) used in v1.2. Because the grid also gives α-element abundances as a free parameter via \[α/Fe\], the metallicity coordinate \[Fe/H\] is no longer exactly equivalent to \[M/H\]. We provide [`BolometricCorrections.MIST.MISTChemistryv2`](@ref) to access the v2.5 solar abundances following the [chemical mixture API](@ref chemistry_api).

```@docs
BolometricCorrections.MIST.MISTChemistryv2
```

The full range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | ≈1,500 K | ≈5,000,000 K |
| logg   | −1.0   | 9.5   |
| \[Fe/H\] | −3.0 dex   | +0.5 dex  |
| \[α/Fe\] | −0.2 dex   | +0.6 dex  |
| Av     | 0.0 mag    | 6.0 mag   |
| Rv     | 3.1    | 3.1   |

!!! note "Teff grid"
    The raw v2.5 tables store ``\log_{10}(T_{\rm eff})`` rather than linear ``T_{\rm eff}``.
    All public API still accepts and returns temperatures in Kelvin.

The full grid of unique values for the dependent variables is available in `BolometricCorrections.MIST.gridinfov2`.

```@example
import BolometricCorrections # hide
keys(BolometricCorrections.MIST.gridinfov2)
```

## Types

```@docs
BolometricCorrections.MIST.MISTBCGridv2
```

The constructor for [`MISTBCGridv2`](@ref) accepts the same human-readable photometric system names as `MISTBCGridv1`. In addition to all photometric systems available in v1.2, v2.5 adds:

- **Euclid** (VIS + NISP)
- **JWST/NIRISS**
- **HST/ACS-SBC** (far-UV channel)
- **RoboAO** (laser guide star AO system)
- **Roman** (Nancy Grace Roman Space Telescope)

!!! note "Roman"
    MIST v1.2 provided BCs for preliminary Roman filters from when it was still called WFIRST. The MIST v2.5 BCs should be preferred for Roman.

The full list of registered photometric systems can be inspected at runtime:

```@example
import DataDeps # hide
using BolometricCorrections.MIST # hide
# custom `show` call to prevent truncation of output # hide
show(stdout, "text/plain", filter(x -> occursin("MIST", x) & occursin("2.5", x), keys(DataDeps.registry)))
```

Once a grid is constructed for a particular photometric system, a BC table with fixed \[Fe/H\], \[α/Fe\], and ``A_V`` can be interpolated.

```@docs
BolometricCorrections.MIST.MISTBCTablev2
```

## Photometric Zeropoints

The MIST v2.5 bolometric corrections adopt the same solar bolometric magnitude ``M_{\odot,\text{bol}} = 4.74`` and bolometric luminosity ``L_\odot = 3.828 \times 10^{33}`` erg s⁻¹ as v1.2. The zeropoint table for converting between the AB, Vega, and ST systems is the same [`BolometricCorrections.MIST.zpt`](@ref) instance documented on the [MIST overview page](@ref MIST).

## MIST v2.5 References
This page cites the following references:

```@bibliography
Pages = ["mist_v2.md"]
Canonical = false
```
