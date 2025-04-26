# [API Reference](@id api)

## BC Grid API
```@docs
BolometricCorrections.AbstractBCGrid
filternames(::BolometricCorrections.AbstractBCGrid)
extrema(::BolometricCorrections.AbstractBCGrid)
Table(::BolometricCorrections.AbstractBCGrid)
columnnames(::BolometricCorrections.AbstractBCGrid)
columns(::BolometricCorrections.AbstractBCGrid)
getproperties(::BolometricCorrections.AbstractBCGrid, ::Tuple{Vararg{Symbol}})
```

## BC Table API
```@docs
BolometricCorrections.AbstractBCTable
filternames(::BolometricCorrections.AbstractBCTable)
extrema(::BolometricCorrections.AbstractBCTable)
Table(::BolometricCorrections.AbstractBCTable)
columnnames(::BolometricCorrections.AbstractBCTable)
columns(::BolometricCorrections.AbstractBCTable)
getproperties(::BolometricCorrections.AbstractBCTable, ::Tuple{Vararg{Symbol}})
```

## [Photometric Zeropoints API](@id zpt_api)
```@docs
BolometricCorrections.AbstractZeropoints
filternames(::BolometricCorrections.AbstractZeropoints)
vegamags
abmags
stmags
BolometricCorrections.Mbol
BolometricCorrections.Lbol
```

## [Chemical Mixture API](@id chemistry_api)
Chemical mixtures used for BCs typically start by defining an assumed chemical mixture for the Sun. When considering the solar chemical mixture, it is important to note that there is a difference between photospheric abundances and protostellar abundances as metals diffuse out of the photosphere over time, leading the photospheric abundances to be lower than protostellar abundances. We differentiate between these in our API. Generally the protostellar quantities (i.e., [`X`](@ref), [`Y`](@ref), [`Z`](@ref)) are preferred for most use cases, as these are the initial conditions for the stellar models set by the researchers and are therefore static. The diffusive processes responsible for the difference between the photospheric and protostellar metallicities depend on the initial mass of a star in addition to its age, and so a population of stars that share a single protostellar metallicity will not all have identical photospheric metallicities, even if they all have the same age. Ideally the photospheric metallicity should be a quantity tracked during the simulation of a stellar model so that it would be included in the stellar track -- the photospheric metallicities are therefore a *prediction* of the stellar evolution models, rather than an intrinsic property set at the beginning of the simulation.

```@docs
BolometricCorrections.AbstractChemicalMixture
chemistry(::BolometricCorrections.AbstractBCTable)
X(::BolometricCorrections.AbstractChemicalMixture)
X_phot(::BolometricCorrections.AbstractChemicalMixture)
Y_p(::BolometricCorrections.AbstractChemicalMixture)
Y(::BolometricCorrections.AbstractChemicalMixture)
Y_phot(::BolometricCorrections.AbstractChemicalMixture)
Z(::BolometricCorrections.AbstractChemicalMixture)
Z_phot(::BolometricCorrections.AbstractChemicalMixture)
```

The above methods define the solar standard assumed in chemical mixture models. The information contained in an `AbstractChemicalMixture` can be used to convert between different chemical conventions (i.e., the metal mass fraction [`Z`](@ref) and logarithmic metallicity [\[M/H\]](@ref MH)). These conversion methods are documented below with accompanying notes for use. 

It has been shown that helium abundance ``Y`` typically only affects broadband BCs at the level of a few thousandths of magnitudes [Girardi2007,VandenBerg2022](@cite). However, it is common for stellar BCs to assume helium abundances that are a function of the stellar metallicity to match the relation used for the stellar evolution models where the helium abundance makes a more significant difference. For grids in which ``Y`` scales with ``Z``, the following methods can be used to derive other quantities.

```@docs
X(::BolometricCorrections.AbstractChemicalMixture, ::Any)
Y(::BolometricCorrections.AbstractChemicalMixture, ::Any)
Z(::BolometricCorrections.AbstractChemicalMixture, ::Any)
MH(::BolometricCorrections.AbstractChemicalMixture, ::Any)
```

Note that we do not offer methods scaling the photospheric abundance values with ``Z``, such as `MH_phot(mix::BolometricCorrections.AbstractChemicalMixture, Z)`, as the diffusive processes that change the photospheric abundances relative to the initial abundances depend on both the age and initial mass of the star in question -- it is therefore inappropriate to simply scale the assumed solar photospheric abundances to other bulk metallicities ``Z``.