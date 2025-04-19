# BolometricCorrections.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgarling.github.io/BolometricCorrections.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgarling.github.io/BolometricCorrections.jl/dev/)
[![Build Status](https://github.com/cgarling/BolometricCorrections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgarling/BolometricCorrections.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgarling/BolometricCorrections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgarling/BolometricCorrections.jl)


**THIS PACKAGE IS PRE-RELEASE**

Provides access to and interpolation of pre-computed grids of bolometric corrections. Currently supported BC grids are

 - [MIST](https://waps.cfa.harvard.edu/MIST/)

Aug 23, 2024: *Rather than have one package (module) that provides different BC libraries as submodules, I think it is better to have separate packages that each provide one BC library. These can share a common interface package like `BolometricCorrectionsBase.jl` or something if necessary, and `BolometricCorrections.jl` could be superpackage that puts all the others together through submodules or something potentially, but for now I'm moving to work on the separate packages for individual libraries.*

Feb 25, 2025: Currently happy with the one package design. There are trade-offs to both single package and multi-package architectures but for now the single package design seems to be sufficient.