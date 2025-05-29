# BolometricCorrections.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgarling.github.io/BolometricCorrections.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgarling.github.io/BolometricCorrections.jl/dev/)
[![Build Status](https://github.com/cgarling/BolometricCorrections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgarling/BolometricCorrections.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgarling/BolometricCorrections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgarling/BolometricCorrections.jl)

A bolometric correction (BC) is the offset between a star's absolute bolometric magnitude and its absolute magnitude in a specific bandpass or filter (e.g., V band). In order to place theoretical stellar models into observational filter spaces, bolometric corrections must be applied to the bolometric magnitudes of the stellar models. Here we provide access to and interpolation of pre-computed grids of [bolometric corrections](https://en.wikipedia.org/wiki/Bolometric_correction). See our documentation linked in the badges above for additional information. Currently supported bolometric correction grids are

 - [MIST](https://waps.cfa.harvard.edu/MIST/)

This package integrates with [StellarTracks.jl](https://github.com/cgarling/StellarTracks.jl) to interpolate isochrones from stellar tracks and apply bolometric corrections to place the isochrones in the observational magnitude space.