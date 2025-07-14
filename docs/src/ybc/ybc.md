# [YBC](@id YBC)

This submodule enables interaction with the "YBC" bolometric correction (BC) grid described in [Chen2019](@citet). These BCs are commonly used in conjunction with the PARSEC stellar models and are used in the "CMD webform" maintained by Leo Girardi. The YBC grid presents a uniform processing of a number of different stellar atmosphere libraries (see section 3.1 of [Chen2019](@citealt)) that each cover different parameter spaces (e.g., surface gravity and effective temperature) as well as different types of stars (e.g., Wolf-Rayet stars, white dwarfs). This approach allows YBC to be robust across a broad range of stellar evolutionary states and makes it very useful for generating synthetic stellar photometry. However, many of these libraries are only used when particular stellar conditions are met, and not all libraries cover the same range of stellar population parameters (metallicity in particular).

Our implementation currently only supports a subset of the libraries in the full YBC set. In particular, we support the [PHOENIX](@ref YBCPHOENIX) models which are used for cool stars, and are working to support the ATLAS9 models used for hotter stars.

## Obtaining Data

Our naming convention for photometric filter systems follows the directory structure of the YBC git repository -- see [supported systems](@ref ybc_systems) below. This repository can be accessed at the URL below.

```@example ybc
using BolometricCorrections.YBC # hide
println(split(YBC.ybc_url, ".git")[1]) # hide
```

The data are obtained by creating a local sparse clone of the YBC git repository with [Scratch.jl](https://github.com/JuliaPackaging/Scratch.jl) and only pulling the data for photometric systems that you request. Presently we only support the standard BCs hosted under the "YBC" subdirectory. These data will be removed automatically if you uninstall the package. In the event you wish to uninstall the data for a particular photometric filter system, you can use [`BolometricCorrections.YBC.remove_table`](@ref), though this should not typically be necessary as the data for each system is on average ~20 MB.

```@docs
BolometricCorrections.YBC.remove_table
```

## [Supported Systems](@id ybc_systems)
Supported filter systems available at the time these docs were built are listed below.

```@setup ybc
function print_in_columns(v::Vector{String}, ncols::Int)
    nrows = ceil(Int, length(v) / ncols)
    padded = vcat(v, fill("", nrows * ncols - length(v)))
    mat = reshape(padded, nrows, ncols)

    # Compute max width per column
    colwidths = [maximum(length.(mat[:, col])) for col in 1:ncols]

    for row in eachrow(mat)
        line = join([rpad(row[i], colwidths[i]) for i in 1:ncols], "  ")  # Adjust spacing here
        println(line)
    end
end
```

```@example ybc
print_in_columns(YBC.systems, 6) # hide
```

## [Chemistry API](@id ybc_chemistry)
We provide the [`BolometricCorrections.YBC.PARSECChemistry`](@ref) type to access information on the solar chemical abundances assumed for the PARSEC models (see also [Bressan2012](@citet)). These solar abundances are also used by the YBC bolometric correction library.

```@docs
BolometricCorrections.YBC.PARSECChemistry
```

Note that in our conversions between ``Z`` and \[M/H\], remembering that `MH = log10(Z/X) - log10(Z⊙/X⊙)`, we use the *photospheric* solar values for `Z⊙` and `X⊙` (these are `Z_⊙` and `X_⊙ = 1 - Z_⊙ - Y_⊙` in Table 3 of [Bressan2012](@citet)). This reproduces the relation between `Z` and \[M/H\] defined in Table 4 of [Bressan2012](@citet), which is also used in the "CMD" webform provided by the PARSEC team.

## YBC References
This page cites the following references:

```@bibliography
Pages = ["ybc.md"]
Canonical = false
```