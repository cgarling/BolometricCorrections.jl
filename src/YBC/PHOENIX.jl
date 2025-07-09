
"""
YBC submodule exposing the bolometric corrections based on PHOENIX BT-Settl atmospheres computed by the YBC team. The original atmosphere models are hosted by SVO [here](https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=bt-settl) and StSCI provides a subset of the PHOENIX library for use with their Synphot software [here](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot). These models assume [Asplund2009](@citet) solar chemical abundances.

Reference articles for these models are [Allard2014](@citet), [Allard2012](@citet)
"""
module PHOENIX

using ..YBC: pull_table, parse_filterinfo
using ArgCheck: @argcheck
using Compat: @compat
using Interpolations: interpolate, extrapolate, Flat, Throw, BSpline, Cubic, Line, OnGrid
# import CSV
using FITSIO: FITS, read_header, colnames
using Printf: @sprintf # Formatted conversion of floats to strings
using TypedTables: Table

""" `NTuple{5, Symbol}` listing the dependent variables in the YBC.PHOENIX BC grid. """
const _dependents = (:logTeff, :logg, :MH, :Av, :Rv)
""" A_v values in files. For each filter "J", each fits file will have columns "J", "J_Av0.5", "J_Av1", and so on."""
const _Av = (0, 0.5, 1, 2, 5, 10, 20)
# const _Av = ["0.5", "1", "2", "5", "10", "20"]
""" Unique values of [M/H] for the PHOENIX BT-Settl models. """
const _mh = (-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.0, 0.3, 0.5)
const _logTeff = range(3.41f0, 4.85f0; step=0.01f0)
const _logg = range(-0.5f0, 6.0f0; step=0.5f0)

# User-facing information on the grid 
""" Unique values for dependent variables in the YBC.PHOENIX bolometric correction grid. """
const gridinfo = (logTeff = _logTeff,
                  logg = _logg,
                  MH = _mh,
                  Av = _Av,
                  Rv = (3.1,))
@compat public gridinfo

"""
    _parse_filename(f::AbstractString)
Return [M/H] and [α/Fe] of PHOENIX BT-Settl model given a filename (example: "Avodonnell94Rv3.1BT-Settl_M-0.0_a+0.0.BC.fits").
"""
function _parse_filename(f::AbstractString)
    f = basename(f)
    @argcheck occursin("BT-Settl", f)
    mh = split(split(f, "_M")[2], "_a")[1]
    α_fe = split(split(f, "_a")[2], ".BC")[1]
    return (MH = parse(Float64, mh), α_fe = parse(Float64, α_fe))
end

function PHOENIXYBCGrid(grid::AbstractString; prefix::AbstractString="YBC")
    path = pull_table(String(grid), String(prefix))
    files = filter(x->occursin("BT-Settl", x), readdir(joinpath(path, "regrid"); join=true))
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names
    return files
    # Dependent variables "logTeff", "logg" should be the same for all
    # BT-Settl models. Therefore to interpolate as a function of [M/H] and Av,
    # we just need to form an array with all the relevant BCs that we can then
    # interpolate between. 
end

function PHOENIXYBCTable(grid::AbstractString, mh::Real, Av::Real; prefix::AbstractString="YBC")
    @argcheck mapreduce(isapprox(mh), |, gridinfo.MH) "Provided [M/H] $mh not in available values $(gridinfo.MH); use YBCPHOENIXGrid for grid interpolation."
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use YBCPHOENIXGrid for grid interpolation."
    path = pull_table(String(grid), String(prefix))
    files = filter(x->occursin("BT-Settl", x), readdir(joinpath(path, "regrid"); join=true))
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names
    # return files

    # Figure out which file we need for given mh
    goodfile = files[findfirst(≈(mh), _parse_filename(file).MH for file in files)]
    # Figure out FITS file column name for given Av
    Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
    # Access FITS file
    FITS(goodfile, "r") do f
        data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
        return data
    end
    # Dependent variables "logTeff", "logg" should be the same for all
    # BT-Settl models. Therefore to interpolate as a function of [M/H] and Av,
    # we just need to form an array with all the relevant BCs that we can then
    # interpolate between. 
end

end # module