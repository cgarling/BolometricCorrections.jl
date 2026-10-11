using BolometricCorrections

# Set up for plotting
import PyPlot as plt
import PyPlot: @L_str # For LatexStrings
plt.rc("text", usetex=true)
plt.rc("font", family="serif", serif=["Computer Modern"], size=16)
plt.rc("figure", figsize=(6,6))
plt.rc("patch", linewidth=1, edgecolor="k", force_edgecolor=true)
# https://matplotlib.org/stable/gallery/images_contours_and_fields/interpolation_methods.html
plt.rc("image", interpolation="none")

function plot_mist_bc_table(gridname::AbstractString, filtername::AbstractString, feh::Number, Av::Number)
    # Load BC grid
    grid = MISTBCGrid(gridname)

    # Create table for fixed [Fe/H], Av
    table = grid(feh, Av)
    table_ext = extrema(table) # keys(table) == (:Teff, :logg)

    # Evaluate BC table on a grid of Teff, logg
    # We'll plot Teff decreasing from left to right
    teff = reverse(range(table_ext.Teff[1], 10_000; length=1000))
    logg = range(table_ext.logg...; length=1000)
    # Transposing logg will evaluate BC for every combination of teff and logg
    data = table.(teff, logg')
    # data will be a (length(teff), length(logg)) matrix, with each element
    # a vector giving the BCs in all the filters. We will pick just one filter
    filter_index = findfirst(==(filtername), String(i) for i in filternames(grid))
    plot_data = [d[filter_index] for d in data]

    fig,ax1 = plt.subplots()
    fig.suptitle("MIST Bolometric Corrections for "*filtername)
    im1 = ax1.imshow(permutedims(plot_data), origin="lower", 
                     extent=(reverse(extrema(teff))..., extrema(logg)...),
                     aspect="auto", rasterized=false)
    ax1.text(0.7, 0.9, "[Fe/H]: "*string(feh)*"\n"*L"A$_V$: "*string(Av)*" mag", transform=ax1.transAxes)

    ax1.set_xlabel(L"T$_{eff}$ [K]")
    ax1.set_ylabel(L"log($g$)")
    fig.colorbar(im1)
    return fig
end
# plot_mist_bc_table("JWST", "F090W", -1, 0)
