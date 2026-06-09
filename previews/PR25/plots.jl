using CairoMakie # CairoMakie re-exports Makie
using BolometricCorrections
using BolometricCorrections: AbstractBCGrid, AbstractBCTable, filternames
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid
set_theme!(theme_latexfonts(); 
           fontsize = 25,
           size = (700, 600),
           Axis = (xticks = Makie.LinearTicks(5),
                   yticks = Makie.LinearTicks(7),
                   xminorticksvisible=true),
            Scatter = (strokecolor=:black, strokewidth=1),
            Colorbar = (ticks = Makie.LinearTicks(7),))
                   # xminorgridvisible=true, 
                   # yminorgridvisible=true))

# For AbstractBCGrid argument, create table and call AbstractBCTable method
# function plot_bc_table(grid::AbstractBCGrid, filtername::AbstractString, mh::Number, Av::Number, args...)
#     table = grid(mh, Av)
#     return plot_bc_table(table, filtername, args...)
# end
# # KoesterWDGrid does not take MH argument
# function plot_bc_table(grid::KoesterWDYBCGrid, filtername::AbstractString, Av::Number, args...)
#     table = grid(Av)
#     return plot_bc_table(table, filtername, args...)
# end
# function plot_bc_table(grid::KoesterWDYBCGrid, filtername::AbstractString, mh::Number, Av::Number, args...)
#     return plot_bc_table(grid, filtername, Av, args...)
# end

# AbstractBCTable method
function plot_bc_table(table::AbstractBCTable, filtername::AbstractString, Teff::AbstractVector, logg::AbstractVector)
    # We'll plot Teff decreasing from left to right
    if !issorted(Teff; rev=true)
        Teff = sort(Teff; rev=true)
    end
    # Transposing logg will evaluate BC for every combination of Teff and logg
    data = table.(Teff, logg')
    logTeff = log10.(Teff)
    filter_index = findfirst(==(filtername), String(i) for i in filternames(table))
    plot_data = [d[filter_index] for d in data]

    f = Figure()
    ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"$\text{log} \ g$")
    ax.xreversed = true
    p = heatmap!(ax, logTeff, logg, plot_data; interpolate=false, colormap=:gist_rainbow) # :viridis) # , colormap=:cividis)
    Colorbar(f[:, end+1], p)

    # # Create top axis to display T [K] x ticks
    # # This messes up the heatmap, don't know why ...
    # # xticks = ax.xaxis.tickvalues[]
    # ax_top = Axis(f[1, 1]; xaxisposition = :top, xlabel = "Teff [K]", xticklabelsvisible = true, xticksvisible = true, yticksvisible=false, yticklabelsvisible=false,
    #               xtickformat = x -> @. string(round(exp10(x), digits=1))) # , xticks = xticks)
    # linkxaxes!(ax, ax_top)
    # # xlims!(ax, reverse(extrema(logTeff)))
    return f, ax
end

# Plots the full Teff, logg extent of the table
function plot_bc_table(table::AbstractBCTable, filtername::AbstractString)
    Teff = logrange(reverse(extrema(table).Teff)...; length=1000)
    logg = range(extrema(table).logg...; length=1000)
    plot_bc_table(table::AbstractBCTable, filtername::AbstractString, Teff, logg)
end

function plot_bc_table_diff(table1::AbstractBCTable, table2::AbstractBCTable, fnames::NTuple{2, String}, Teff::AbstractVector, logg::AbstractVector)
    # We'll plot Teff decreasing from left to right
    if !issorted(Teff; rev=true)
        Teff = sort(Teff; rev=true)
    end
    logTeff = log10.(Teff)
    # Transposing logg will evaluate BC for every combination of Teff and logg
    data1 = table1.(Teff, logg')
    filter_index1 = findfirst(==(fnames[1]), String(i) for i in filternames(table1))
    plot_data1 = [d[filter_index1] for d in data1]
    # Repeat for second table
    data2 = table2.(Teff, logg')
    filter_index2 = findfirst(==(fnames[2]), String(i) for i in filternames(table2))
    plot_data2 = [d[filter_index2] for d in data2]

    f = Figure()
    ax = Axis(f[1, 1], xlabel="log(Teff)", ylabel=L"$\text{log} \ g$")
    ax.xreversed = true
    p = heatmap!(ax, logTeff, logg, plot_data1 .- plot_data2; interpolate=false, colormap=:gist_rainbow) # :viridis) # , colormap=:cividis)
    Colorbar(f[:, end+1], p)
    return f, ax
end