# BolometricCorrections.jl
Interfaces for bolometric corrections libraries.

*Rather than have one package (module) that provides different BC libraries as submodules, I think it is better to have separate packages that each provide one BC library. These can share a common interface package like `BolometricCorrectionsBase.jl` or something if necessary, and `BolometricCorrections.jl` could be superpackage that puts all the others together through submodules or something potentially, but for now I'm moving to work on the separate packages for individual libraries.*