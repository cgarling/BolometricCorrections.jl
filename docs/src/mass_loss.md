# [Stellar Mass Loss Rates](@id mass_loss)

Some BC grids for high-mass stars include the stellar mass-loss (outflow) rate as an additional dependent parameter (e.g., the **YBC WM-basic grid**). To use these grids, we must be able to calculate stellar mass loss rates from a star's other parameters as mass-loss rates are not typically tracked self-consistently in stellar interior models. Most implementations employ scaling relations between the mass-loss rate and parameters like the stellar metallicity and luminosity for radiation-driven winds typical in high-mass O and B stars. For example, PARSEC uses multiple different scaling relations to determine the stellar mass-loss rate depending on the phase of stellar evolution, including [Vink2001](@citet) -- see section 2.3 of [Tang2014](@citet) for more details. We include a few such models below.

## Björklund 2021
```@docs
BolometricCorrections.Bjorklund2021MassLoss
```

Below we reproduce Figure 7 from [Bjorklund2021](@citet). 
```@example plotting
using CairoMakie
using BolometricCorrections: Bjorklund2021MassLoss
model = Bjorklund2021MassLoss()
Zvals = [model.Zsol, 0.5 * model.Zsol, 0.2 * model.Zsol]
labels = ["Solar", "LMC", "SMC"]
logL = range(6-1.4, 6; length=100)
Mdot = model.(Zvals, logL') # Makes matrix with size (length(Zvals), length(logL))

fig = Figure()
# ax = Axis(fig[1, 1], xlabel=L"\log\left(\frac{L}{10^6 \, L_\odot}\right)", ylabel=L"\log \left( \frac{\dot{\text{M}}}{\text{M}_\odot \, \text{yr}^{-1}} \right)", yticks = -9:1:-6, xticks = -1.4:0.2:0.0) # hide
ax = Axis(fig[1, 1], xlabel=L"\log (\text{L} \ \left[ \text{L}_\odot \right]) - 6", ylabel=L"\log \left( \dot{\text{M}} \ \left[ \text{M}_\odot \ \text{yr}^{-1} \right] \right)", yticks = -9:1:-6, xticks = -1.4:0.2:0.0)
# Label(fig[2, 1], L"\log \left(\frac{L}{10^6 \, L_\odot}\right)"; halign = :center, valign = :bottom, padding = (0, 0, 5, 0), tellwidth = false) # hide
for i in eachindex(Zvals, labels)
    lines!(ax, logL .- 6, log10.(Mdot[i, :]), label=labels[i])
end
axislegend(ax, position=:lt)
fig
```


## Mass Loss References
This page cites the following references:

```@bibliography
Pages = ["mass_loss.md"]
Canonical = false
```