
"""
YBC submodule exposing the bolometric corrections based on PHOENIX BT-Settl atmospheres computed by the YBC team. The original atmosphere models are hosted by SVO [here](https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=bt-settl) and StSCI provides a subset of the PHOENIX library for use with their Synphot software [here](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot). These models assume [Asplund2009](@citet) solar chemical abundances.

The main reference article for these models is [Allard2012](@citet). [Allard2013](@citet) discusses BT-Settl models with the solar composition from [Caffau2011](@citet), but this module uses the model with the [Asplund2009](@citet) abundances.
"""
module PHOENIX

using ...BolometricCorrections: repack_submatrix, AbstractBCTable, AbstractBCGrid, interp1d, interp2d
import ...BolometricCorrections: zeropoints, filternames, chemistry # Y_p, X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry, vegamags, abmags, stmags, Mbol, Lbol
using ...BolometricCorrections.MIST: MISTChemistry # MIST and YBC PHOENIX both use Asplund2009 abundances, so just use MISTChemistry
using ..YBC: pull_table, parse_filterinfo, check_prefix

using ArgCheck: @argcheck
using Compat: @compat
# using Interpolations: interpolate, extrapolate, Flat, Throw, BSpline, Cubic, Line, OnGrid
using Interpolations: cubic_spline_interpolation, Throw, Flat
# import CSV
using FITSIO: FITS, read_header, colnames
using Printf: @sprintf # Formatted conversion of floats to strings
using TypedTables: Table
using StaticArrays: SVector

export PHOENIXYBCTable, PHOENIXYBCGrid, chemistry, filternames

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
""" Unique values for dependent variables in the YBC.PHOENIX bolometric correction grid. Note that `logg = 6` is missing for [M/H] = -2.5, -3, -3.5, -4 """
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
# AbstractBCGrid code


struct PHOENIXYBCGrid{A <: Number, C <: AbstractVector{A}, N} <: AbstractBCGrid{A}
    data::Matrix{Matrix{A}} # A should be Float32
    mag_zpt::C
    systems::Vector{String}
    name::String
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end

function PHOENIXYBCGrid(data::Matrix{Matrix{A}}, mag_zpt::AbstractArray{<:Number}, systems, name, filternames) where {A}
    return PHOENIXYBCGrid(data, convert.(A, mag_zpt), String.(systems), String(name), tuple(Symbol.(filternames)...))
end

function PHOENIXYBCGrid(grid::AbstractString; prefix::AbstractString="YBC")
    check_prefix(prefix)
    path = pull_table(String(grid), String(prefix))
    files = filter(x->occursin("BT-Settl", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    elseif length(files) != length(gridinfo.MH)
        error("Number of files found for grid $grid is $(length(files)), expected $(length(gridinfo.MH)). Data may be corrupted. \\
        Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names
    # mh_vals = [_parse_filename(file).MH for file in files]
    # Sort files by [M/H] value
    idxs = sortperm([_parse_filename(file).MH for file in files])
    files = files[idxs]
    # Check that MH values as parsed from filenames match gridinfo.MH
    if ~mapreduce(isapprox, &, gridinfo.MH, _parse_filename(file).MH for file in files)
        error("File [M/H] values not as expected -- please report.")
    end
    # Dependent variables "logTeff", "logg" should be the same for all
    # BT-Settl models. Therefore to interpolate as a function of [M/H] and Av,
    # we just need to form an array with all the relevant BCs that we can then
    # interpolate between.
    # Read all data and pack into dense matrix
    data = Matrix{Matrix{Float32}}(undef, length(files), length(gridinfo.Av))
    for i in eachindex(files)
        file = files[i]
        FITS(file, "r") do f
            for j in eachindex(gridinfo.Av)
                Av = gridinfo.Av[j]
                # Figure out FITS file column name for given Av
                Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(Av)
                data[i,j] = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
            end
        end
    end
    return PHOENIXYBCGrid(data, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, filternames)
end
(grid::PHOENIXYBCGrid)(mh::Real, Av::Real) = PHOENIXYBCTable(grid, mh, Av)
Base.show(io::IO, z::PHOENIXYBCGrid) = print(io, "YBC PHOENIX bolometric correction grid for photometric system $(z.name).")
# function Table(grid::PHOENIXYBCGrid)
#     data = grid.data
#     tables = Vector{Table}(undef, length(data))
# end
Base.extrema(::PHOENIXYBCGrid) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                  logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                  MH = (first(gridinfo.MH), last(gridinfo.MH)),
                                  Av = (first(gridinfo.Av), last(gridinfo.Av)),
                                  Rv = (first(gridinfo.Rv), last(gridinfo.Rv)))
filternames(grid::PHOENIXYBCGrid) = grid.filters
chemistry(::Type{<:PHOENIXYBCGrid}) = MISTChemistry()
# zeropoints(::MISTBCGrid) = zpt


#########################################################
# A single BC table, with fixed [M/H] and Av
struct PHOENIXYBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    MH::A
    Av::A
    mag_zpt::Vector{A}
    systems::Vector{String}
    name::String
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
function PHOENIXYBCTable(MH::Real, Av::Real, mag_zpt::Vector{<:Real}, systems, name, itp, filters)
    T = promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
    return PHOENIXYBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), itp, filters)
end
chemistry(::Type{<:PHOENIXYBCTable}) = MISTChemistry()
Base.show(io::IO, z::PHOENIXYBCTable) = print(io, "YBC PHOENIX BT-Settl bolometric correction table with for system $(z.name) with [M/H] ",
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
    grid, prefix = String(grid), String(prefix)
    check_prefix(prefix)
    @argcheck mapreduce(isapprox(mh), |, gridinfo.MH) "Provided [M/H] $mh not in available values $(gridinfo.MH); use YBCPHOENIXGrid for grid interpolation."
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use YBCPHOENIXGrid for grid interpolation."
    path = pull_table(grid, prefix)
    files = filter(x->occursin("BT-Settl", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    elseif length(files) != length(gridinfo.MH)
        error("Number of files found for grid $grid is $(length(files)), expected $(length(gridinfo.MH)). Data may be corrupted. \\
        Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Figure out which file we need for given [M/H]
    goodfile = files[findfirst(≈(mh), _parse_filename(file).MH for file in files)]
    # Figure out FITS file column name for given Av
    Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
    # Access FITS file
    FITS(goodfile, "r") do f
        data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames) # ← 900μs ↑ 354.021 μs ↓ 160 μs
        # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
        # Account for the fact that logg = 6 is missing for mh <= -2.5
        println(length(unique(read(f[2], "logg"))) * length(unique(read(f[2], "logTeff"))))
        return data
        logg = mh > -2.5 ? gridinfo.logg : gridinfo.logg[begin:end-1]
        newdata = repack_submatrix(data, length(logg), length(gridinfo.logTeff), Val(length(filternames)))
        itp = cubic_spline_interpolation((logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
        return PHOENIXYBCTable(mh, Av, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
    end
end

function PHOENIXYBCTable(grid::PHOENIXYBCGrid, mh::Real, Av::Real)
    check_vals(mh, Av)
    filters = filternames(grid)
    data = grid.data

    Av_vec = SVector(gridinfo.Av) # Need vectors to use searchsortedfirst
    MH_vec = SVector(gridinfo.MH)

    # Exact values are in grid; no interpolation necessary
    if mh ∈ gridinfo.MH && Av ∈ gridinfo.Av
        # Account for the fact that logg = 6 is missing for mh <= -2.5.
        logg = mh > -2.5 ? gridinfo.logg : gridinfo.logg[begin:end-1]
        submatrix = data[searchsortedfirst(MH_vec, mh), searchsortedfirst(Av_vec, Av)]
    else
        if mh ∈ gridinfo.MH
            # Account for the fact that logg = 6 is missing for mh <= -2.5.
            logg = mh > -2.5 ? gridinfo.logg : gridinfo.logg[begin:end-1]
            MH_idx = searchsortedfirst(MH_vec, mh)
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            mat1 = data[MH_idx, Av_idx]
            mat2 = data[MH_idx, Av_idx + 1]
            submatrix = interp1d(Av, Av_vec[Av_idx], Av_vec[Av_idx + 1], mat1, mat2)
        elseif Av ∈ gridinfo.Av
            # Account for the fact that logg = 6 is missing for mh <= -2.5.
            logg = mh > -2 ? gridinfo.logg : gridinfo.logg[begin:end-1]
            Av_idx = searchsortedfirst(Av_vec, Av)
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mat1 = data[MH_idx, Av_idx]
            mat2 = data[MH_idx + 1, Av_idx]
            # Any interpolation between -2.5 < mh < -2 
            # will have data matrices of different shapes -- properly truncate longer matrix
            if -2.5 < mh < -2
                # logg is iterated first in mat1 and mat2; dims are (logg, logTeff). 
                # Need to skip every Nth element, where N = length(gridinfo.logg)
                filtered = [i for i in axes(mat2, 1) if i % length(gridinfo.logg) != 0]
                mat2 = mat2[filtered, :]
            end
            submatrix = interp1d(mh, MH_vec[MH_idx], MH_vec[MH_idx + 1], mat1, mat2)
        else
            # Account for the fact that logg = 6 is missing for mh <= -2.5.
            logg = mh > -2 ? gridinfo.logg : gridinfo.logg[begin:end-1]
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            Av1, Av2 = Av_vec[Av_idx], Av_vec[Av_idx + 1]
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mh1, mh2 = MH_vec[MH_idx], MH_vec[MH_idx + 1]

            mat1_1 = data[MH_idx, Av_idx]
            mat2_1 = data[MH_idx + 1, Av_idx]
            mat1_2 = data[MH_idx, Av_idx + 1]
            mat2_2 = data[MH_idx + 1, Av_idx + 1]
            if -2.5 < mh < -2 # See above
                filtered = [i for i in axes(mat2_1, 1) if i % length(gridinfo.logg) != 0]
                mat2_1 = mat2_1[filtered, :]
                mat2_2 = mat2_2[filtered, :]
            end
            # Perform bilinear interpolation
            submatrix = interp2d(mh, Av, mh1, mh2, Av1, Av2, mat1_1, mat2_1, mat1_2, mat2_2)
        end
    end
    newdata = repack_submatrix(submatrix, length(logg), length(gridinfo.logTeff), Val(length(filters)))
    itp = cubic_spline_interpolation((logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
    return PHOENIXYBCTable(mh, Av, grid.mag_zpt, grid.systems, grid.name, itp, filters)
end

end # module