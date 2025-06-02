# BolometricCorrections.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgarling.github.io/BolometricCorrections.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgarling.github.io/BolometricCorrections.jl/dev/)
[![Build Status](https://github.com/cgarling/BolometricCorrections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgarling/BolometricCorrections.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgarling/BolometricCorrections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgarling/BolometricCorrections.jl)

A bolometric correction (BC) is the offset between a star's absolute bolometric magnitude and its absolute magnitude in a specific bandpass or filter (e.g., V band). In order to place theoretical stellar models into observational filter spaces, bolometric corrections must be applied to the bolometric magnitudes of the stellar models. Here we provide access to and interpolation of pre-computed grids of [bolometric corrections](https://en.wikipedia.org/wiki/Bolometric_correction). See our documentation linked in the badges above for additional information. Currently supported bolometric correction grids are

 - [MIST](https://waps.cfa.harvard.edu/MIST/)

This package integrates with [StellarTracks.jl](https://github.com/cgarling/StellarTracks.jl) to interpolate isochrones from stellar tracks and apply bolometric corrections to place the isochrones in the observational magnitude space.

## Installation

This package is registered to Julia's General registry and can be installed via Pkg from the Julia REPL by executing

```julia
import Pkg;
Pkg.add("BolometricCorrections");
```

## Example Usage

The HST ACS WFC bolometric corrections from the MIST grid can be loaded with

```julia
using BolometricCorrections
grid = MISTBCGrid("hst_acs_wfc")
```

If you haven't yet acquired these data, you will be prompted to allow the data to be downloaded and installed.

The MIST bolometric correction grid covers a wide range of metallicity and reddening (`Av`) values. We can interpolate the grid at specific values of metallicity ([M/H]) and reddening (`Av`), obtaining a `MISTBCTable` with

```julia
mh, Av = -1.01, 0.125
table = grid(mh, Av)
```

This table can now be interpolated at any supported value of effective temperature (`Teff`, in Kelvin) and surface gravity (`logg`) as

```julia
Teff, logg = 2755, 0.01
table(Teff, logg)
```

which will return a vector containing the bolometric corrections in all the filters defined in the HST ACS WFC system. These filters can be listed with

```julia
filternames(grid)
```
```
(:ACS_WFC_F435W, :ACS_WFC_F475W, :ACS_WFC_F502N, :ACS_WFC_F550M, :ACS_WFC_F555W, :ACS_WFC_F606W, :ACS_WFC_F625W, :ACS_WFC_F658N, :ACS_WFC_F660N, :ACS_WFC_F775W, :ACS_WFC_F814W, :ACS_WFC_F850LP, :ACS_WFC_F892N)
```

To enable convenient calculation of BCs over many `Teff`, `logg` entries, you can call `table(Teff, logg)` with array arguments,

```julia
table([2755, 2756], [0.01, 0.02])
```

which return the BCs stacked into a matrix. You can also get the results returned as a `TypedTables.Table` with

```julia
using TypedTables: Table
bcs = table(Table, [2755, 2756], [0.01, 0.02])
```
```
Table with 13 columns and 2 rows:
     ACS_WFC_F435W  ACS_WFC_F475W  ACS_WFC_F502N  ACS_WFC_F550M  ACS_WFC_F555W  ACS_WFC_F606W  ACS_WFC_F625W  ACS_WFC_F658N  ⋯
   ┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 1 │ -6.79721       -5.9089        -6.67736       -4.90698       -5.21007       -4.28292       -3.98601       -3.32683       ⋯
 2 │ -6.77482       -5.88998       -6.65945       -4.88979       -5.19314       -4.2692        -3.973         -3.31317       ⋯
 ```

so that the results can be more conveniently accessed with the following syntax

```julia
bcs.ACS_WFC_F435W # Retrieve BCs in F435W filter
```
