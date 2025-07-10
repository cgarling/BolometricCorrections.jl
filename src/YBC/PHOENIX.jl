
"""
YBC submodule exposing the bolometric corrections based on PHOENIX BT-Settl atmospheres computed by the YBC team. The original atmosphere models are hosted by SVO [here](https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=bt-settl) and StSCI provides a subset of the PHOENIX library for use with their Synphot software [here](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot). These models assume [Asplund2009](@citet) solar chemical abundances.

The main reference article for these models is [Allard2012](@citet). [Allard2013](@citet) discusses BT-Settl models with the solar composition from [Caffau2011](@citet), but this module uses the model with the [Asplund2009](@citet) abundances.
"""
module PHOENIX

using ...BolometricCorrections: repack_submatrix, AbstractBCTable
import ...BolometricCorrections: zeropoints, filternames, chemistry # Y_p, X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry, vegamags, abmags, stmags, Mbol, Lbol
using ...BolometricCorrections.MIST: MISTChemistry # MIST and YBC PHOENIX both use Asplund2009 abundances, so just use MISTChemistry
using ..YBC: pull_table, parse_filterinfo, check_prefix

using ArgCheck: @argcheck
using Compat: @compat
# using Interpolations: interpolate, extrapolate, Flat, Throw, BSpline, Cubic, Line, OnGrid
using Interpolations: cubic_spline_interpolation, Throw
# import CSV
using FITSIO: FITS, read_header, colnames
using Printf: @sprintf # Formatted conversion of floats to strings
using TypedTables: Table

export PHOENIXYBCTable

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

"""
    check_vals(mh, Av)

Validate that [M/H] value `mh` and ``A_V`` value `Av` are valid for the YBC PHOENIX BC grid.
Throws `DomainError` if check fails, returns `nothing` if check is successful.
"""
function check_vals(mh, Av)
    mh_ext = extrema(gridinfo.MH)
    if mh < first(mh_ext) || mh > last(mh_ext)
        throw(DomainError(mh, "Provided [M/H] $mh is outside the bounds for the YBC PHOENIX BC tables $mh_ext"))
    end
    Av_ext = extrema(gridinfo.Av)
    if Av < first(Av_ext) || Av > last(Av_ext)
        throw(DomainError(Av, "Provided A_v $Av is outside the bounds for the YBC PHOENIX BC tables $Av_ext"))
    end
end

#########################################################
# A single BC table, with fixed [M/H] and Av
struct PHOENIXYBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    MH::A
    Av::A
    mag_zpt::Vector{A}
    system::Vector{String}
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
function PHOENIXYBCTable(MH::Real, Av::Real, mag_zpt::Vector{<:Real}, system, itp, filters)
    T = promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
    return PHOENIXYBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, system), itp, filters)
end
chemistry(::PHOENIXYBCTable) = MISTChemistry()
Base.show(io::IO, z::PHOENIXYBCTable) = print(io, "YBC PHOENIX BT-Settl bolometric correction table with [M/H] ",
                                          z.MH, " and V-band extinction ", z.Av)
filternames(table::PHOENIXYBCTable) = table.filters
zeropoints(table::PHOENIXYBCTable) = table.mag_zpt
# Interpolations uses `bounds` to return interpolation domain
# We will just use the hard-coded grid bounds; extremely fast
Base.extrema(::PHOENIXYBCTable) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                   logg = (first(gridinfo.logg), last(gridinfo.logg)))
# Base.extrema(::PHOENIXYBCTable) = (Teff = extrema(exp10.(gridinfo.logTeff)), logg = extrema(gridinfo.logg))
(table::PHOENIXYBCTable)(Teff::Real, logg::Real) = table.itp(logg, log10(Teff))
# Data are naturally Float32 -- convert Float64 args for faster evaluation (~35% faster)
(table::PHOENIXYBCTable)(Teff::Float64, logg::Float64) = table.itp(convert(Float32, logg), log10(convert(Float32, Teff)))
# to broadcast over both teff and logg, you do table.(teff, logg')

function PHOENIXYBCTable(grid::AbstractString, mh::Real, Av::Real; prefix::AbstractString="YBC")
    check_prefix(prefix)
    @argcheck mapreduce(isapprox(mh), |, gridinfo.MH) "Provided [M/H] $mh not in available values $(gridinfo.MH); use YBCPHOENIXGrid for grid interpolation."
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use YBCPHOENIXGrid for grid interpolation."
    path = pull_table(String(grid), String(prefix))
    files = filter(x->occursin("BT-Settl", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Figure out which file we need for given [M/H]
    goodfile = files[findfirst(≈(mh), _parse_filename(file).MH for file in files)]
    # Figure out FITS file column name for given Av
    Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
    # Access FITS file
    FITS(goodfile, "r") do f
        data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
        # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
        newdata = repack_submatrix(data, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filternames)))
        itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Throw())
        return PHOENIXYBCTable(mh, Av, filterinfo.mag_zeropoint, string.(filterinfo.photometric_system), itp, tuple(Symbol.(filternames)...))
    end
end

# function PHOENIXYBCGrid(grid::AbstractString; prefix::AbstractString="YBC")
#     check_prefix(prefix)
#     path = pull_table(String(grid), String(prefix))
#     files = filter(x->occursin("BT-Settl", x), readdir(joinpath(path, "regrid"); join=true))
#     if length(files) == 0
#         error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
#     end
#     filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
#     filternames = filterinfo.names
#     return files
#     # Dependent variables "logTeff", "logg" should be the same for all
#     # BT-Settl models. Therefore to interpolate as a function of [M/H] and Av,
#     # we just need to form an array with all the relevant BCs that we can then
#     # interpolate between.
# end

end # module