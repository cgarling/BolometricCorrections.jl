using CairoMakie # CairoMakie re-exports Makie
using BolometricCorrections: AbstractBCGrid, AbstractBCTable, filternames
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid
set_theme!(theme_latexfonts(); 
           fontsize = 20,
           size = (600, 600),
           Axis = (xticks = Makie.LinearTicks(5),
                   yticks = Makie.LinearTicks(7),
                   # xticks=Makie.WilkinsonTicks(10; k_min=5, k_max=5),
                   # yticks=Makie.WilkinsonTicks(5; k_min=5, k_max=5)))
                   # xminorticks=Makie.IntervalsBetween(5),
                   xminorticksvisible=true),
            Scatter = (strokecolor=:black, strokewidth=1))
                   # xminorgridvisible=true, 
                   # yminorgridvisible=true))

# For AbstractBCGrid argument, create table and call AbstractBCTable method
function plot_bc_table(grid::AbstractBCGrid, filtername::AbstractString, mh::Number, Av::Number, args...)
    table = grid(mh, Av)
    return plot_bc_table(table, filtername, args...)
end
# KoesterWDGrid does not take MH argument
function plot_bc_table(grid::KoesterWDYBCGrid, filtername::AbstractString, Av::Number, args...)
    table = grid(Av)
    return plot_bc_table(table, filtername, args...)
end
function plot_bc_table(grid::KoesterWDYBCGrid, filtername::AbstractString, mh::Number, Av::Number, args...)
    return plot_bc_table(grid, filtername, Av, args...)
end

# AbstractBCTable method
function plot_bc_table(table::AbstractBCTable, filtername::AbstractString, Teff::AbstractVector, logg::AbstractVector)
    # We'll plot Teff decreasing from left to right
    println("running")
    if !issorted(Teff; rev=true)
        Teff = sort(Teff; rev=true)
    end
    # Transposing logg will evaluate BC for every combination of Teff and logg
    data = table.(Teff, logg')
    filter_index = findfirst(==(filtername), String(i) for i in filternames(table))
    plot_data = [d[filter_index] for d in data]

    f = Figure()
    ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"$\text{log} \ g$")
    ax.xreversed = true
    p = heatmap!(ax, log10.(Teff), logg, plot_data; interpolate=false, colormap=:viridis) # , colormap=:cividis)
    Colorbar(f[:, end+1], p)
    return f, ax
end

# Plots the full Teff, logg extent of the table
function plot_bc_table(table::AbstractBCTable, filtername::AbstractString)
    Teff = logrange(reverse(extrema(table).Teff)...; length=1000)
    logg = range(extrema(table).logg...; length=1000)
    plot_bc_table(table::AbstractBCTable, filtername::AbstractString, Teff, logg)
end