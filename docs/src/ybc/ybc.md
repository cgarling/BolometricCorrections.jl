# [YBC](@id YBC)

This submodule enables interaction with the "YBC" bolometric correction (BC) grid described in [Chen2019](@citet). These BCs are commonly used in conjunction with the PARSEC stellar models and are used in the "CMD webform" maintained by Leo Girardi. The YBC grid presents a uniform processing of a number of different stellar atmosphere libraries (see section 3.1 of [Chen2019](@citealt)) that each cover different parameter spaces (e.g., surface gravity and effective temperature) as well as different types of stars (e.g., Wolf-Rayet stars, white dwarfs). This approach allows YBC to be robust across a broad range of stellar evolutionary states and makes it very useful for generating synthetic stellar photometry. However, many of these libraries are only used when particular stellar conditions are met, and not all libraries cover the same range of stellar population parameters (metallicity in particular).

Our implementation currently only supports a subset of the libraries in the full YBC set. In particular, we support the [PHOENIX](@ref YBCPHOENIX) models (used for cool stars), the [ATLAS9](@ref YBCATLAS9) models (used for hot stars), the [Koester & Tremblay white dwarf models](@ref YBCKoester), and the [WM-basic O and B star models](@ref YBCWMbasic). Interpolation between these libraries as a function of stellar properties is supported by the [`YBCGrid`](@ref) and [`YBCTable`](@ref) types.

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

We illustrate which model library is used as a function of surface gravity `logg` and effective temperature `Teff` for stars without outflows (`Mdot = 0`) below.

```@example ybc
using BolometricCorrections # hide
using BolometricCorrections.YBC # hide
include(joinpath(@__DIR__, "..", "plots.jl")) # hide
grid = YBCGrid("jwst_nircam_wide") # hide
# Teff = range(extrema(grid).Teff[1], 10_000; length=1000) # hide
Teff = logrange(exp10(3.5), 10_000; length=1000) # hide
logg = range(extrema(grid).logg...; length=1000) # hide
f, ax = plot_bc_table(grid(-1, 0), "F090W", Teff, logg) # hide
ax.title = "YBC BCs for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text=L"[M/H] = -1\n Av = 0\n$\dot{M}=0$", align=(:right, :top), space=:relative) # hide
vlines!(ax, log10.([grid.transitions.koester.Teff[2],]), color = "black") # hide
hlines!(ax, [grid.transitions.koester.logg[1],], color = "black"; # hide
        xmin=(log10(grid.transitions.koester.Teff[2]) - log10(minimum(Teff)))/diff(log10.([minimum(Teff), maximum(Teff)]))[1]) # hide
# lines!(ax, repeat([log10(grid.transitions.koester.Teff[2])], 2), # hide
#        [grid.transitions.koester.logg[2], maximum(logg)], color=:black, linestyle=:dash) # hide
# lines!(ax, log10.([grid.transitions.koester.Teff[2], maximum(Teff)]), # hide
#        repeat([grid.transitions.koester.logg[2]], 2), color=:black, linestyle=:dash) # hide
text!(ax, 0.05, 0.95, text="Koester\nWhite dwarfs", align=(:left, :top), space=:relative) # hide
text!(ax, 0.05, 0.35, text="ATLAS9\nHot stars", align=(:left, :top), space=:relative) # hide
text!(ax, 0.95, 0.35, text="PHOENIX\nCool stars", align=(:right, :top), space=:relative) # hide
f # hide
```

We illustrate the transition between the ATLAS9 models used for hot stars and the WM-basic models used for O- and B-type stars with radiation-driven outflows at fixed surface gravity below.

```@example ybc
Teff = logrange(18_000, 24_000; length=1000) # hide
logTeff = log10.(Teff) # hide
logg = 2.5 # hide
Mdot = logrange(1e-10, 1e-4; length=1000) # hide
logMdot = log10.(Mdot) # hide

table = grid(-1, 0) # hide
data = table.(Teff, logg, Mdot') # hide

filter_index = findfirst(==("F090W"), String(i) for i in filternames(table)) # hide
plot_data = [d[filter_index] for d in data] # hide

f = Figure() # hide
ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"\text{log} \left( \dot{M} \right)") # hide
ax.xreversed = true # hide
p = heatmap!(ax, logTeff, logMdot, plot_data; interpolate=false, colormap=:gist_rainbow) # :viridis) # , colormap=:cividis) # hide
Colorbar(f[:, end+1], p) # hide
text!(ax, 0.95, 0.95, text=L"[M/H] = -1\n Av = 0\n$\text{log} \, g=%$logg$", align=(:right, :top), space=:relative) # hide
vlines!(ax, log10.([grid.transitions.wmbasic.Teff[1],]), color = "black", # hide
        ymin=(log10(grid.transitions.wmbasic.Mdot[2]) - log10(minimum(Mdot)))/diff(log10.([minimum(Mdot), maximum(Mdot)]))[1]) # hide
hlines!(ax, [log10(grid.transitions.wmbasic.Mdot[2]),], color = "black"; # hide
        xmin=(log10(grid.transitions.wmbasic.Teff[1]) - log10(minimum(Teff)))/diff(log10.([minimum(Teff), maximum(Teff)]))[1]) # hide
text!(ax, 0.05, 0.95, text="WM-basic", align=(:left, :top), space=:relative) # hide
text!(ax, 0.05, 0.35, text="ATLAS9", align=(:left, :top), space=:relative) # hide

f # hide
```

## Types

```@docs
YBCGrid
YBCTable
```

## [Chemistry API](@id ybc_chemistry)
We provide the [`BolometricCorrections.YBC.PARSECChemistry`](@ref) type to access information on the solar chemical abundances assumed for the PARSEC stellar models (see also [Bressan2012](@citet)) and YBC bolometric correction library.

```@docs
BolometricCorrections.YBC.PARSECChemistry
```

Note that in our conversions between ``Z`` and \[M/H\], remembering that `MH = log10(Z/X) - log10(Z⊙/X⊙)`, we use the *photospheric* solar values for `Z⊙` and `X⊙` (these are `Z_⊙` and `X_⊙ = 1 - Z_⊙ - Y_⊙` in Table 3 of [Bressan2012](@citet)). This reproduces the relation between `Z` and \[M/H\] defined in Table 4 of [Bressan2012](@citet), which is also used in the "CMD" webform provided by the PARSEC team.

The individual submodules that constitute the YBC library have different assumptions about the solar chemical mixture. When you interpolate a [`YBCGrid`](@ref) to a particular metallicity \[M/H\], this is converted to the corresponding metal mass fraction [`Z`](@ref) for [`PARSECChemistry`](@ref BolometricCorrections.YBC.PARSECChemistry). Each of the submodules are then interpolated to this common metal mass fraction -- in this way, all the submodules are normalized to the same metal mass fraction, though their solar chemical mixtures are not identical. This is the same approach taken by [Chen2019](@citet).

## YBC References
This page cites the following references:

```@bibliography
Pages = ["ybc.md"]
Canonical = false
```