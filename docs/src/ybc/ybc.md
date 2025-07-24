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
ybcgrid = YBCGrid("jwst_nircam_wide") # hide
# Teff = range(extrema(ybcgrid).Teff[1], 10_000; length=1000) # hide
Teff = logrange(exp10(3.5), 10_000; length=1000) # hide
logg = range(extrema(ybcgrid).logg...; length=1000) # hide
f, ax = plot_bc_table(ybcgrid(-1, 0), "F090W", Teff, logg) # hide
ax.title = "YBC BCs for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text=L"[M/H] = -1\n Av = 0\n$\dot{M}=0$", align=(:right, :top), space=:relative) # hide
vlines!(ax, log10.([ybcgrid.transitions.koester.Teff[2],]), color = "black") # hide
hlines!(ax, [ybcgrid.transitions.koester.logg[1],], color = "black"; # hide
        xmin=(log10(ybcgrid.transitions.koester.Teff[2]) - log10(minimum(Teff)))/diff(log10.([minimum(Teff), maximum(Teff)]))[1]) # hide
# lines!(ax, repeat([log10(ybcgrid.transitions.koester.Teff[2])], 2), # hide
#        [ybcgrid.transitions.koester.logg[2], maximum(logg)], color=:black, linestyle=:dash) # hide
# lines!(ax, log10.([ybcgrid.transitions.koester.Teff[2], maximum(Teff)]), # hide
#        repeat([ybcgrid.transitions.koester.logg[2]], 2), color=:black, linestyle=:dash) # hide
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

table = ybcgrid(-1, 0) # hide
data = table.(Teff, logg, Mdot') # hide

filter_index = findfirst(==("F090W"), String(i) for i in filternames(table)) # hide
plot_data = [d[filter_index] for d in data] # hide

f = Figure() # hide
ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"\text{log} \left( \dot{M} \right)") # hide
ax.xreversed = true # hide
p = heatmap!(ax, logTeff, logMdot, plot_data; interpolate=false, colormap=:gist_rainbow) # :viridis) # , colormap=:cividis) # hide
Colorbar(f[:, end+1], p) # hide
text!(ax, 0.95, 0.95, text=L"[M/H] = -1\n Av = 0\n$\text{log} \, g=%$logg$", align=(:right, :top), space=:relative) # hide
vlines!(ax, log10.([ybcgrid.transitions.wmbasic.Teff[1],]), color = "black", # hide
        ymin=(log10(ybcgrid.transitions.wmbasic.Mdot[2]) - log10(minimum(Mdot)))/diff(log10.([minimum(Mdot), maximum(Mdot)]))[1]) # hide
hlines!(ax, [log10(ybcgrid.transitions.wmbasic.Mdot[2]),], color = "black"; # hide
        xmin=(log10(ybcgrid.transitions.wmbasic.Teff[1]) - log10(minimum(Teff)))/diff(log10.([minimum(Teff), maximum(Teff)]))[1]) # hide
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

## [Metallicity Extrapolation](@id ybc_extrapolation)
The sub-libraries that make up the YBC grid do not have uniform metallicity coverage. The metallicity coverage for the PHOENIX, ATLAS9, and WM-basic libraries is shown below. The Koester & Tremblay white dwarf library does not consider metallicity as it assumes pure hydrogen atmospheres for the white dwarfs.

```@example ybc
NamedTuple{(:phoenix, :atlas9, :wmbasic)}(extrema(grid).MH for grid in (PHOENIXYBCGrid, ATLAS9YBCGrid, BolometricCorrections.YBC.WMbasic.WMbasicYBCGrid))
```

As such, there is a choice to be made regarding how to treat metallicities that fall outside the range of one or more of the constituent libraries. We generally wish to maximize the useful extent of the BC grid while minimizing errors due to extrapolation. For cool stars, PHOENIX is used which has the greatest metallicity range. Attempts to interpolate a [`YBCTable`](@ref) outside the bounds of the PHOENIX metallicity range will throw an error. Extrapolation behavior is controlled by the `extrapolate::Bool` keyword argument that can be passed to [`YBCGrid`](@ref) which is `true` by default. When `extrapolate = true`, the ATLAS9 and WM-basic models will be extrapolated with a flat boundary condition to the same metallicity range as supported by the PHOENIX libary. Details and justification for this behavior is given below. This extrapolation can be disabled by setting the keyword `extrapolate = false` when constructing a [`YBCGrid`](@ref).

### ATLAS9

For hot stars, we can examine the differences in the BCs at the edges of the metallicity grids to determine whether the BCs are changing rapidly or not (essentially estimating the rate of change of the BCs with respect to the metallicity). The differences in the ATLAS9 BCs for the JWST/NIRCam F090W filter between metallicities \[M/H\] = -2.5 and \[M/H\] = -2.4 are shown below. The dynamic range of the differences is quite low, indicating that the BCs are not changing rapidly as the metallicity decreases -- the magnitude of the derivative of the bolometric correction with respect to the metallicity is at most d(BC) / d(\[M/H\]) = 0.024 mag / dex while the median is closer to 0.01 mag / dex. 

```@example ybc
Teff = logrange(6_000, 50_000; length=1000) # hide
logTeff = log10.(Teff) # hide
logg = range(extrema(ATLAS9YBCGrid).logg...; length=1000) # hide

atlas9grid = ATLAS9YBCGrid("jwst_nircam_wide") # hide
table1, table2 = atlas9grid(-2.5, 0.0), atlas9grid(-2.4, 0.0) # hide
# table1, table2 = ATLAS9YBCGrid("jwst_nircam_wide")(-2.5, 0.0), PHOENIXYBCGrid("jwst_nircam_wide")(-2.5, 0.0) # hide

f, ax = plot_bc_table_diff(table1, table2, ("F090W", "F090W"), Teff, logg) # hide
ax.title = "ATLAS9 F090W [M/H]=$(MH(table1)) minus [M/H]=$(MH(table2))" # hide
text!(ax, 0.95, 0.95, text="Av = 0", align=(:right, :top), space=:relative) # hide

f # hide
```

By way of comparison, the difference between the ATLAS9 BCs and the PHOENIX BCs at \[M/H\] = -2.5 can be much larger:

```@example ybc
Teff = logrange(6_000, 50_000; length=1000) # hide
logTeff = log10.(Teff) # hide
logg = range(extrema(ATLAS9YBCGrid).logg...; length=1000) # hide

phoenixgrid = PHOENIXYBCGrid("jwst_nircam_wide") # hide
table1, table2 = atlas9grid(-2.5, 0.0), phoenixgrid(-2.5, 0.0) # hide

f, ax = plot_bc_table_diff(table1, table2, ("F090W", "F090W"), Teff, logg) # hide
ax.title = L"ATLAS9 $-$ PHOENIX, F090W, [M/H]=%$(MH(table1))" # hide
text!(ax, 0.95, 0.95, text="Av = 0", align=(:right, :top), space=:relative) # hide

f # hide
```

As the ATLAS9 models generally show better agreement with observations for hot stars, our default implementation (`extrapolate = true` in the [`YBCGrid`](@ref) constructor) is to extrapolate the ATLAS9 bolometric grid to match the metallicity range of the PHOENIX grid. This extrapolation is "flat" -- when a metallicity outside the range of the ATLAS9 grid is requested (e.g., \[M/H\] = -3), we will use the BC table for the nearest valid metallicity (\[M/H\] = -2.5). Note that the upper limit of the metallicity grid for PHOENIX and ATLAS9 are the same (\[M/H\] = 0.5) so no extrapolation is used at the high-metallicity end.

### WM-basic

We can perform the same examination for the WM-basic models, which have a more limited range of metallicity. Here we examine the change between the WM-basic BCs for the JWST/NIRCam F090W filter with fixed ``\dot{M} = 10^{-6}`` solar masses per year (the mean value, as the grid covers ``10^{-7}`` to ``10^{-5}``). The median rate of change of the BC with the metallicity is still fairly low d(BC) / d(\[M/H\]) = 0.05 mag / dex, but clearly there are some regions of parameter space where larger variation is exhibited.

```@example ybc
using BolometricCorrections.YBC.WMbasic: WMbasicYBCGrid # hide
Teff = logrange(6_000, 50_000; length=1000) # hide
logTeff = log10.(Teff) # hide
logg = range(extrema(WMbasicYBCGrid).logg...; length=1000) # hide
Mdot = 1e-6 # hide
mh_min = minimum(extrema(WMbasicYBCGrid).MH) # hide

wmbasic_grid = WMbasicYBCGrid("jwst_nircam_wide")
table1, table2 = wmbasic_grid(mh_min, 0.0), wmbasic_grid(mh_min + 0.1, 0.0) # hide

# Transposing logg will evaluate BC for every combination of Teff and logg # hide
data1 = table1.(Teff, logg', Mdot) # hide
filter_index1 = findfirst(==("F090W"), String(i) for i in filternames(table1)) # hide
plot_data1 = [d[filter_index1] for d in data1] # hide
# Repeat for second table # hide
data2 = table2.(Teff, logg', Mdot) # hide
filter_index2 = findfirst(==("F090W"), String(i) for i in filternames(table2)) # hide
plot_data2 = [d[filter_index2] for d in data2] # hide

f = Figure() # hide
ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"$\text{log} \ g$") # hide
ax.xreversed = true # hide
p = heatmap!(ax, logTeff, logg, plot_data1 .- plot_data2; interpolate=false, colormap=:gist_rainbow) # :viridis) # , colormap=:cividis) # hide
Colorbar(f[:, end+1], p) # hide

ax.title = "WM-basic, F090W, [M/H]=$(round(MH(table1), digits=2)) minus [M/H]=$(round(MH(table2), digits=2))" # hide
text!(ax, 0.95, 0.95, text=L"Av = 0\n$\text{log} \left( \dot{M} \right) = %$(log10(Mdot))$", align=(:right, :top), space=:relative) # hide

f # hide
```

Given that the WM-basic models are for hot O- and B-type stars, the alternative in YBC would be to use the ATLAS9 models which we discussed above. The ATLAS9 models have wider metallicity coverage, but do not consider the effects of outflowing material, which leads to significant discrepancies.

```@example ybc
Teff = logrange(6_000, 50_000; length=1000) # hide
logTeff = log10.(Teff) # hide
logg = range(extrema(ATLAS9YBCGrid).logg...; length=1000) # hide
Mdot = 1e-6 # hide
mh_min = minimum(extrema(WMbasicYBCGrid).MH) # hide

table1, table2 = wmbasic_grid(mh_min, 0.0), atlas9grid(mh_min, 0.0) # hide

# Transposing logg will evaluate BC for every combination of Teff and logg # hide
data1 = table1.(Teff, logg', Mdot) # hide
filter_index1 = findfirst(==("F090W"), String(i) for i in filternames(table1)) # hide
plot_data1 = [d[filter_index1] for d in data1] # hide
# Repeat for second table # hide
data2 = table2.(Teff, logg') # hide
filter_index2 = findfirst(==("F090W"), String(i) for i in filternames(table2)) # hide
plot_data2 = [d[filter_index2] for d in data2] # hide

f = Figure() # hide
ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"$\text{log} \ g$") # hide
ax.xreversed = true # hide
p = heatmap!(ax, logTeff, logg, plot_data1 .- plot_data2; interpolate=false, colormap=:gist_rainbow) # :viridis) # , colormap=:cividis) # hide
Colorbar(f[:, end+1], p) # hide

ax.title = L"WM-basic $-$ ATLAS9, F090W, [M/H]=%$(round(MH(table1), digits=2))" # hide
text!(ax, 0.95, 0.95, text=L"Av = 0\n$\text{log} \left( \dot{M} \right) = %$(log10(Mdot))$", align=(:right, :top), space=:relative) # hide

f # hide
```

Clearly error due to extrapolating the WM-basic models in metallicity will be significantly less than the error we would otherwise incur by using the ATLAS9 models, which do not include the effects of stellar winds. Therefore, when using [`YBCGrid`](@ref) with keyword argument `extrapolate = true`, we extrapolate the WM-basic grid with flat boundary conditions to match the metallicity coverage of the PHOENIX library.


## YBC References
This page cites the following references:

```@bibliography
Pages = ["ybc.md"]
Canonical = false
```