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
Chemical mixtures used for BCs typically start by defining an assumed chemical mixture for the Sun. When considering the solar chemical mixture, it is important to note that there is a difference between photospheric abundances and protostellar abundances as metals diffuse out of the photosphere over time, leading the photospheric abundances to be lower than protostellar abundances. We differentiate between these in our API. Generally the protostellar quantities (i.e., [`X`](@ref), [`Y`](@ref), [`Z`](@ref)) are preferred for most use cases.

```@docs
BolometricCorrections.AbstractChemicalMixture
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
X_phot(::BolometricCorrections.AbstractChemicalMixture, ::Any)
Y(::BolometricCorrections.AbstractChemicalMixture, ::Any)
Y_phot(::BolometricCorrections.AbstractChemicalMixture, ::Any)
Z(::BolometricCorrections.AbstractChemicalMixture, ::Any)
Z_phot(::BolometricCorrections.AbstractChemicalMixture, ::Any)
MH(::BolometricCorrections.AbstractChemicalMixture, ::Any)
```