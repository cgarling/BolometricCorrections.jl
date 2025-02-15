module MIST

# using ..BolometricCorrections: Table, columnnames # relative path for parent module
using ..BolometricCorrections: AbstractBCGrid, AbstractBCTable
import ..BolometricCorrections: filternames
using CodecXz: XzDecompressorStream # Decompress downloaded BCs
import CSV
using DataDeps: register, DataDep, @datadep_str
# using DataInterpolations: PCHIPInterpolation # AbstractInterpolation, AkimaInterpolation, LinearInterpolation, CubicSpline, CubicHermiteSpline,
using Interpolations: interpolate, Linear, Gridded # scale, BSpline
# import HDF5 # not currently using
using Printf: @sprintf # Formatted conversion of floats to strings
using StaticArrays: SVector
import Tar
import Tables # for Tables.matrix conversion
import TypedTables: Table, columnnames, getproperties, columns

# Initialization and datadeps
include("init.jl")

""" Unique values of `Teff` in the MIST BC tables. """
const _mist_Teff = (2500.0, 2800.0, 3000.0, 3200.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, 60000.0, 70000.0, 80000.0, 90000.0, 100000.0, 110000.0, 120000.0, 130000.0, 140000.0, 150000.0, 160000.0, 170000.0, 180000.0, 190000.0, 200000.0, 300000.0, 400000.0, 500000.0, 600000.0, 700000.0, 800000.0, 900000.0, 1.0e6) 
""" Unique values of `logg` in the MIST BC tables. """
const _mist_logg = (-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5)
""" Unique values of `Av` in the MIST BC tables. """
const _mist_Av = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0)
""" Unique values of `Rv` in the MIST BC tables. """
const _mist_Rv = (3.1,)
""" Unique values of [Fe/H] in the MIST BC tables. """
const _mist_feh = (-4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75)

function check_vals(feh, Av)
    feh_ext = extrema(_mist_feh)
    if feh < first(feh_ext) || feh > last(feh_ext)
        throw(DomainError(feh, "Provided [Fe/H] $feh is outside the bounds for the MIST BC tables $feh_ext"))
    end
    Av_ext = extrema(_mist_Av)
    if Av < first(Av_ext) || Av > last(Av_ext)
        throw(DomainError(Av, "Provided A_v $Av is outside the bounds for the MIST BC tables $Av_ext"))

    end
end

# Read a MIST BC file post-processed in init.jl into Table
function read_mist_bc_processed(fname::AbstractString)
    return CSV.read(fname, Table)
end
# Given a datadep path, return path to post-processed file
mist_processed_fname(fname::AbstractString) = joinpath(fname, last(splitpath(fname))*".gz")

# Data structures
# A is default data type -- inferred from provided table
struct MISTBCGrid{A,B,C <: AbstractString} <: AbstractBCGrid{A}
    table::B # Usually a TypedTables.Table
    filename::C
    # unique_feh::D
    # unique_Av::D
    # unique_Rv::D # MIST only has 1 Rv
    function MISTBCGrid(table::B, filename::C) where {B, C}
        A = Base.promote_eltype(first(table))
        # D = Vector{A}
        new{A, B, C}(table, filename)
    end
end
function MISTBCGrid(grid::AbstractString)
    if mapreduce(x->occursin(x,grid), |, ("JWST", "jwst"))
        fname = mist_processed_fname(datadep"JWST")
        return MISTBCGrid(read_mist_bc_processed(fname), fname)
    end
end
Table(grid::MISTBCGrid) = grid.table
# columnnames(grid::MISTBCGrid) = columnnames(Table(grid))
# columns(grid::MISTBCGrid) = columns(Table(grid))
# getproperties(grid::MISTBCGrid, names::Tuple{Vararg{Symbol}}) = getproperties(Table(grid), names) 
# A function that will extract the dependent variables from a MIST BC grid
const _mist_dependents = (:Teff, :logg, :feh, :Av, :Rv)
extract_dependents(grid::MISTBCGrid) = getproperties(grid, _mist_dependents)
Base.extrema(grid::MISTBCGrid) = NamedTuple{_mist_dependents}(extrema(col) for col in columns(extract_dependents(grid)))
# filternames(grid::MISTBCGrid) = [string(name) for name in columnnames(grid)[6:end]]
filternames(grid::MISTBCGrid) = columnnames(grid)[length(_mist_dependents)+1:end]


#########################################################
# A single BC table, with fixed feh and Av
struct MISTBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    feh::A
    Av::A
    # Rv::A # Only one unique RV for MIST 
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} # NTuple{N, Symbol} giving filter names
end
filternames(table::MISTBCTable) = table.filters

# MISTBCTable(feh::Real, Av::Real, Rv::Real, itp, filters) = MISTBCTable(promote(feh, Av, Rv)..., itp, filters)
# MISTBCTable(feh::Real, Av::Real, itp, filters) = MISTBCTable(promote(feh, Av)..., itp, filters)

# Use statically known size from filters argument to repack submatrix
# into a vector of SVectors to pass into interpolator
function _repack_submatrix(submatrix::AbstractArray{T},
                           filters::NTuple{N, Symbol}) where {T, N}
    submatrix = reshape(submatrix, length(_mist_logg),
                        length(_mist_Teff),
                        length(filters))
    return [SVector{N, T}(view(submatrix,i,j,:)) for i=axes(submatrix,1),j=axes(submatrix,2)]
end
function MISTBCTable(feh::Real, Av::Real, grid::MISTBCGrid)
    # ext = extrema(grid)
    # @assert ext.feh[1] ≤ feh ≤ ext.feh[2]
    check_vals(feh, Av)

    # Exact values are in grid; no interpolation necessary
    if feh ∈ _mist_feh && Av ∈ _mist_Av
        table = Table(grid)
        # Basically all the time is spent in this filter ... 
        subtable = filter(row -> (row.feh ≈ feh) && (row.Av ≈ Av), table)
        filters = filternames(grid)
        submatrix = Tables.matrix(getproperties(subtable, filters))
        # submatrix = reshape(submatrix, 26, 70, 29)
        # submatrix = reshape(submatrix, length(_mist_logg),
        #                     length(_mist_Teff),
        #                     length(filters))
        # return [SVector{29, Float64}(view(submatrix,i,j,:)) for i=axes(submatrix,1),j=axes(submatrix,2)]

        # return submatrix, (subtable.Teff, subtable.logg)
        # BSpline only supports ranges for x and y, not general ...
        # itp = scale(interpolate(submatrix, BSpline(Linear())), (subtable.Teff, subtable.logg))
        # itp = interpolate((subtable.Teff, subtable.logg),
        # itp = interpolate((unique(subtable.logg), unique(subtable.Teff)),
        #                   [submatrix[i,j,:] for i=axes(submatrix,1),j=axes(submatrix,2)],
        #                   Gridded(Linear()))
        # Make use of constant information
        itp = interpolate((collect(_mist_logg), collect(_mist_Teff)),
                          # [submatrix[i,j,:] for i=axes(submatrix,1),j=axes(submatrix,2)],
                          # [SVector{29, Float64}(view(submatrix,i,j,:)) for i=axes(submatrix,1),j=axes(submatrix,2)], # 2x faster than vector
                          _repack_submatrix(submatrix, filters),
                          Gridded(Linear()))
        
        return MISTBCTable(feh, Av, itp, filters)
        
    end
end
# _name_itp_tuple(val::T, filters::F)::NamedTuple{F,T} where {T,F} = NamedTuple{filters}(val)
# _name_itp_tuple(filters::F, val::T) where {F,T} = NamedTuple{filters}(val)
# _name_itp_tuple(val, filters::Val{T}) where T = NamedTuple{T}(val)
# function (table::MISTBCTable{A,B,N})(Teff::Real, logg::Real) where {A,B,N}
#     # return NamedTuple{table.filters}(Tuple(table.itp(logg, Teff)))
#     # return Tuple(table.itp(logg, Teff))
#     val = table.itp(logg, Teff)
#     filters = filternames(table)
#     # return NamedTuple{filters, NTuple{N,A}}(Tuple(table.itp(logg, Teff)))#::NamedTuple{filters, NTuple{N,A}}
#     return _name_itp_tuple(Tuple(val), Val(filters))
#     # return _name_itp_tuple(filters, Tuple(val))
# end
(table::MISTBCTable)(Teff::Real, logg::Real) = table.itp(logg, Teff)
# to broadcast over both teff and logg, you do table.(teff, logg')

# Not happy to be using Dierckx but you can use Interpolations.jl
# but the results are worse than Dierckx because of the continuity
# boundary conditions; see 
# https://discourse.julialang.org/t/difference-between-two-interpolation-packages-when-evaluating-derivatives-interpolations-jl-and-dierckx-jl/61747


end # module
