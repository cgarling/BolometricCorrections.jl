# [WM-basic](@id YBCWMbasic)

This submodule (`BolometricCorrections.YBC.WMbasic`) enables interaction with the "YBC" bolometric correction (BC) grid of hot star (O and B-type) atmospheres calculated with the WM-basic code by [Chen2015](@citet) -- see section 3.2 in particular. Section 3.3 of [Chen2019](@citet) discusses the calculation of BCs from these models. Radiation-driven outflows can significantly impact the SEDs for stars of these types, and the models are calculated for three different values of mass outflow rate: ``\text{log} (\dot{M}) = (-7, -6, -5)`` where ``\dot{M}`` is in units of solar masses per year. The full grid coverage is given below.

|        | min    | max   |
|--------|--------|-------|
| Teff   | 19,952 K | 100,000 K |
| logg   | 2.5  | 6.0   |
| Mdot   | 1e-7 | 1e-5 |
| Av     | 0.0 mag    | 20.0 mag   |
| Rv     | 3.1    | 3.1   |

!!! warning
    These models **cannot** be used without specifying a mass outflow rate (`Mdot`). Not all stellar evolutionary libraries provide this quantity. Even when they do, these outflow rates are often based on scaling relations with other stellar evolutionary quantities (e.g., metallicity, luminosity) rather than being self-consistently inferred from the interior model. We provide utilities to estimate outflow rates from other stellar quantities using scaling relations [here](@ref mass_loss).

```@example
using BolometricCorrections # hide
using BolometricCorrections.YBC.WMbasic: WMbasicYBCGrid # hide
include(joinpath(@__DIR__, "..", "plots.jl")) # hide
grid = WMbasicYBCGrid("jwst_nircam_wide") # hide
Teff = logrange(exp10(4.25), 100_000; length=1000) # hide
logg = range(extrema(grid).logg...; length=1000) # hide
# f, ax = plot_bc_table(grid(-1, 0), "F090W", Teff, logg) # hide
table = grid(-1, 0) # hide
Mdot = 1e-5 # hide
data = table.(Teff, logg', Mdot) # hide
filter_index = findfirst(==("F090W"), String(i) for i in filternames(table)) # hide
plot_data = [d[filter_index] for d in data] # hide

f = Figure() # hide
ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"$\text{log} \ g$") # hide
ax.xreversed = true # hide
p = heatmap!(ax, log10.(Teff), logg, plot_data; interpolate=false, colormap=:viridis) # , colormap=:cividis) # hide
Colorbar(f[:, end+1], p) # hide

ax.title = "YBC WM-basic BCs for JWST/NIRCam F090W" # hide
text!(ax, 0.95, 0.95, text=L"$A_V = 0$\n$\log\left(\dot{M}\right) = %$(log10(Mdot))$", align=(:right, :top), space=:relative) # hide
f # hide
```

## Types

```@docs
BolometricCorrections.YBC.WMbasic.WMbasicYBCGrid
BolometricCorrections.YBC.WMbasic.WMbasicYBCTable
```

## Chemistry

These models use the same solar chemical abundances as are used for the PARSEC stellar evolutionary library [Bressan2012](@citep), which is represented by the [`BolometricCorrections.YBC.PARSECChemistry`](@ref) type.

## WM-basic References
This page cites the following references:

```@bibliography
Pages = ["wm_basic.md"]
Canonical = false
```