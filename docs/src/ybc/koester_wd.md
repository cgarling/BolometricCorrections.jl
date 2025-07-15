# [Koester & Tremblay White Dwarfs](@id YBCKoester)

This submodule (`BolometricCorrections.YBC.KoesterWD`) enables interaction with the "YBC" bolometric correction (BC) grid calculated from the white dwarf libraries of [Koester2010](@citet) and [Tremblay2009](@citet). These works used Gaia DR2 parallaxes to obtain absolute photometry for white dwarfs, enabling better calibration of the appropriate BCs. The model atmospheres are hosted by SVO [here](http://svo2.cab.inta-csic.es/theory/newov2/index.php). They provide the following description of the models:

!!! info 
    These models are for white dwarfs of spectral type DA with pure hydrogen atmospheres. They use LTE (local thermodynamic equilibrium), hydrostatic equilibrium and plane-parallel, one-dimensional structure. Basic methods and data are described in Koester (2010, Mem.S.A.It. Vol. 81, 921). Since then many improvements were implemented, most notably the hydrogen Stark profiles by Tremblay & Bergeron (2009, ApJ 696,1755), and Tremblay (2015, priv. comm).

Details on the BC calculations are given in section 3.6 of [Chen2019](@citet). As these models should *only* be used for white dwarfs, we do not export the grid and table types from this submodule, indicating that they are treated as internal. This library of models is only used for white dwarf stars in the integrated **YBC grid** implementation, with details on **that page**. These models do not vary with metallicity as they consider only pure hydrogen atmospheres. The full range of dependent variables covered by this grid is given in the table below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 5,011 K | 79,432 K |
| logg   | 6.5  | 9.5   |
| Av     | 0.0 mag    | 20.0 mag   |
| Rv     | 3.1    | 3.1   | 

```@example
using BolometricCorrections # hide
using BolometricCorrections: KoesterWDYBCGrid # hide
include(joinpath(@__DIR__, "..", "plots.jl")) # hide
grid = KoesterWDYBCGrid("jwst_nircam_wide") # hide
Teff = logrange(exp10(3.7), 10_000; length=1000) # hide
logg = range(extrema(grid).logg...; length=1000) # hide
f, ax = plot_bc_table(grid(0), "F090W", Teff, logg) # hide
ax.title = "YBC ATLAS9 BCs for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text="Av = 0", align=(:right, :top), space=:relative) # hide
f # hide
```

## Types

```@docs
BolometricCorrections.YBC.KoesterWD.KoesterWDYBCGrid
BolometricCorrections.YBC.KoesterWD.KoesterWDYBCTable
```

## Chemistry
No chemistry information is provided as these models assume pure hydrogen atmospheres. As such, calling [`chemistry`](@ref) with one of the above types will return `missing`.

```@example
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid, chemistry
chemistry(KoesterWDYBCGrid("acs_wfc"))
```

## KoesterWD References
This page cites the following references:

```@bibliography
Pages = ["koester_wd.md"]
Canonical = false
```