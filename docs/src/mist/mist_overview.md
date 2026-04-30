# [MIST](@id MIST)

This submodule enables interaction with the bolometric correction (BC) grids released as part of the Mesa Isochrones & Stellar Tracks ([MIST](https://mist.science); [Dotter2016,Choi2016,Dotter2026,Bauer2026](@citet)) project. The MIST BC grids are convenient for stellar population work: they cover a wide range of photometric filters, span the full range of effective temperature, surface gravity, and metallicity relevant for stellar evolution, and are *regular* in their dependent variables, which greatly simplifies interpolation.

Two versions of the MIST BC tables are currently supported:

- **[MIST v1.2](@ref MIST_v1)** — the original grid with BCs tabulated as a function of temperature (Teff), surface gravity (logg), metallicity (\[Fe/H\]), and V-band extinction (Av).
- **[MIST v2.5](@ref MIST_v2)** — new grid adding \[α/Fe\] as a free parameter, updated solar chemical abundances, and additional photometric systems.

## Common features

Both versions share:

- The same set of photometric zeropoints stored in [`BolometricCorrections.MIST.zpt`](@ref), an instance of [`BolometricCorrections.MIST.MISTZeropoints`](@ref). This can be used to convert bolometric corrections between the AB, Vega, and ST photometric systems.
- The same solar bolometric magnitude ``M_{\odot,\text{bol}} = 4.74`` and bolometric luminosity ``L_\odot = 3.828 \times 10^{33}`` erg s⁻¹.
- A regular Cartesian grid in all dependent variables, enabling simple bilinear (or higher-order) interpolation.
- The same [call signature API](@ref api): construct a *grid* object by loading a photometric system, then call it to obtain a *table* with the stellar-population parameters fixed.

## Version comparison

| Feature | v1.2 | v2.5 |
|---------|------|------|
| Solar abundances | [Asplund2009](@citet) | [Grevesse1998](@citet) |
| α-element enhancement | fixed (scaled-solar) | \[α/Fe\] ∈ {−0.2, 0, 0.2, 0.4, 0.6} |
| Teff grid | linear (2,500 – 1,000,000 K; 70 pts) | log₁₀ (∼1,500 – 5,000,000 K; 41 pts) |
| \[Fe/H\] range | −4.0 to +0.75 (18 values) | −3.0 to +0.5 (15 values) |
| logg range | −4.0 to 9.5 (26 values) | −1.0 to 9.5 (22 values) |
| Photometric systems | 22 | 28 (adds Euclid, NIRISS, RoboAO, Roman, HST/ACS-SBC) |
| Chemistry type | [`MISTv1Chemistry`](@ref BolometricCorrections.MIST.MISTv1Chemistry) | [`MISTv2Chemistry`](@ref BolometricCorrections.MIST.MISTv2Chemistry) |
| Grid type | [`MISTv1BCGrid`](@ref BolometricCorrections.MIST.MISTv1BCGrid) | [`MISTv2BCGrid`](@ref BolometricCorrections.MIST.MISTv2BCGrid) |
| Table type | [`MISTv1BCTable`](@ref BolometricCorrections.MIST.MISTv1BCTable) | [`MISTv2BCTable`](@ref BolometricCorrections.MIST.MISTv2BCTable) |

Note that for ease of transition from with earlier versions of BolometricCorrections.jl that only supported MIST v1.2, the following deprecated aliases are available

 - `MISTChemistry` aliases for [`MISTv1Chemistry`](@ref BolometricCorrections.MIST.MISTv1Chemistry)
 - `MISTBCTable` aliases for [`MISTv1BCTable`](@ref)
 - `MISTBCGrid` aliases for [`MISTv1BCGrid`](@ref)

## Shared types

The zeropoint table and its type are shared between both versions:

```@docs
BolometricCorrections.MIST.zpt
BolometricCorrections.MIST.MISTZeropoints
```

## MIST References
This page cites the following references:

```@bibliography
Pages = ["mist_overview.md"]
Canonical = false
```
