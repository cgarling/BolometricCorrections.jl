# [Overview](@id overview)

This package enables interaction with grids and tables of stellar [bolometric corrections](https://en.wikipedia.org/wiki/Bolometric_correction). Different grids/libraries are supported through submodules. Currently supported model grids with accompanying submodules are

 - [Mesa Isochrones and Stellar Tracks (MIST)](@ref MIST) [Dotter2016,Choi2016](@cite)
 - [YBC](@ref YBC) [Chen2019](@cite), which interpolates between several different BC libaries -- supported sub-libraries are
   - [PHOENIX](@ref YBCPHOENIX)
   - [ATLAS9](@ref YBCATLAS9)
   - [Koester & Tremblay white dwarf library](@ref YBCKoester)

## [Background](@id background)
A **bolometric correction** is the offset between a star's absolute bolometric magnitude ``M_\text{bol}`` and its absolute magnitude in a specific bandpass or filter ``\lambda``,

```math
BC_{\lambda} = M_\text{bol} - M_\lambda
```

As stars of different types exhibit different stellar spectra, ``BC_{\lambda}`` depends on the properties of the star under consideration. In order to produce mock stellar populations in observational bandpasses, we can calculate these BCs from synthetic stellar spectra -- see section 2.2 of [Girardi2002](@citet) for details on these calculations.

Because the calculation of BCs is expensive, we do not compute new BCs from scratch when we want to produce new mock stellar populations from theoretical stellar models. Instead, large grids of BCs are computed across the relevant stellar parameters once and interpolated to obtain intermediate values as necessary. **These grids of pre-computed BCs are the focus of this package.**

The dependent variables varied across the grid are typically effective temperature ``\left(T_{eff}\right)``, surface gravity ``\left(\text{log}\left(g\right)\right)``, and metallicity (metal mass fraction *Z*, \[M/H\], or \[Fe/H\] depending on convention). Some grids will additionally add interstellar reddening (usually quantified by V-band extinction ``A_V``) and/or offer other chemical compositions beyond solar-scaled chemical mixtures (typically α-enhanced or α-depleted mixtures). Grids of BCs typically differ in their coverage of these dependent variables and the synthetic spectra / stellar atmospheres they use to calculate the BCs. It is therefore useful to have access to several different libraries of BCs so that researchers can quantify how their analyses are affected by the choice of BC library. **These types of investigations are what we hope to facilitate with this package.**

## Photometric Zeropoints
The normalization of the bolometric corrections as defined above relies on the measurement of the solar luminosity. An arbitrary absolute bolometric magnitude ``M_\text{bol}`` is chosen to correspond to this luminosity, setting the normalization. In 2015 the IAU adopted [Resolution B2](@cite Mamajek2015) arguing for a bolometric magnitude normalization defined so that ``M_\text{bol}=0`` corresponds to ``L_\text{bol} = 3.0128 \times 10^{28}`` watts or ``3.0128 \times 10^{35} \, \text{erg} \, \text{s}^{-1}``. With such a normalization, the solar luminosity ``L_\odot = 3.828 \times 10^{26}`` watts corresponds with an absolute bolometric magnitude of ``M_{\text{bol}, \odot} = 4.74``. This normalization has been widely adopted following the IAU's resolution, but older BCs may follow different conventions. The documentation for each BC grid will contain information on the photometric zeropoint assumed to normalize the BCs. 

## Home References
This page cites the following references:

```@bibliography
Pages = ["index.md"]
Canonical = false
```

Also see the [full bibliography](@ref bibliography) for further references cited in this documentation.