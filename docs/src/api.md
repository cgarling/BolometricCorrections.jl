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