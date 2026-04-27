"""
Submodule enabling interaction with the MIST grid of stellar bolometric corrections.
"""
module MIST

# using ..BolometricCorrections: Table, columnnames # relative path for parent module
using ..BolometricCorrections: AbstractBCGrid, AbstractBCTable, AbstractZeropoints, AbstractChemicalMixture,
                               interp1d, interp2d, repack_submatrix
import ..BolometricCorrections: zeropoints, filternames, vegamags, abmags, stmags, Mbol, Lbol, Y_p, X, X_phot,
                                Y, Y_phot, Z, Z_phot, MH, chemistry, gridname

using ArgCheck: @argcheck
using CodecXz: XzDecompressorStream # Decompress downloaded BCs
using Compat: @compat # for @compat public <x>
import CSV
using DataDeps: register, DataDep, @datadep_str
using DataDeps: registry # contains list of registered datadeps
# using DataInterpolations: PCHIPInterpolation # AbstractInterpolation, AkimaInterpolation, LinearInterpolation, CubicSpline, CubicHermiteSpline,
using Interpolations: interpolate, Linear, Gridded, extrapolate, Flat # scale, BSpline
# import HDF5 # not currently using
using Printf: @sprintf # Formatted conversion of floats to strings
using StaticArrays: SVector
import Tar
import Tables # for Tables.matrix conversion
import TypedTables: Table # import to extend
using TypedTables: columnnames, getproperties, columns
using Unicode: normalize # To normalize string arguments

export MISTBCGridv1, MISTBCTablev1, MISTBCGridv2, MISTBCTablev2,
       MISTBCGrid, MISTBCTable,  # deprecated aliases for v1 types
       MISTChemistry, MISTChemistryv1, MISTChemistryv2,
       X, X_phot, Y, Y_phot, Z, Z_phot, Y_p, MH

"`NTuple{5, Symbol}` listing the dependent variables in the MIST v1.2 BC grid."
const _mist_v1_dependents = (:Teff, :logg, :feh, :Av, :Rv)
"`NTuple{6, Symbol}` listing the dependent variables in the MIST v2.5 BC grid (adds `afe` for [α/Fe])."
const _mist_v2_dependents = (:lgTef, :logg, :feh, :afe, :Av, :Rv)
"`NTuple{5, Symbol}` giving iteration order (inner→outer) for v1.2 post-processed tables."
const _mist_v1_dependents_order = (:logg, :Teff, :Av, :Rv, :feh)
"`NTuple{6, Symbol}` giving iteration order (inner→outer) for v2.5 post-processed tables."
const _mist_v2_dependents_order = (:logg, :lgTef, :Av, :afe, :feh, :Rv)

""" Struct to contain the MIST zeropoint information.
A constant instance is available as [`BolometricCorrections.MIST.zpt`](@ref).
Instances are callable with a `filtername::AbstractString` which will return
the zeropoint entry for the relevant filter. The full table can be retrieved
with `Table(zpt)` and the list of available filters can be retrieved with
`filternames(zpt)`. """
struct MISTZeropoints{T} <: AbstractZeropoints
    table::T
end
function Base.getproperty(zpt::MISTZeropoints, name::Symbol)
    # table = Table(zpt)
    if name == :table
        return getfield(zpt, :table)
    end
    # filters = filternames(zpt)
    # strname = String(name)
    # for (i, filt) in enumerate(filters)
    #     if strname == filt
    #         return Table(zpt)[i]
    #     end
    # end
    # println(filters)
    # throw(ArgumentError("""No exact match found for $filtername in MIST filterset, please review above list of filters."""))
    return zpt(String(name))
end
Table(zpt::MISTZeropoints) = zpt.table
Base.show(io::IO, zpt::MISTZeropoints) = display(Table(zpt))
filternames(zpt::MISTZeropoints) = Table(zpt).filter
function (zpt::MISTZeropoints)(filtername::AbstractString)
    filters = filternames(zpt)
    for (i, filt) in enumerate(filters)
        if filtername == filt
            return Table(zpt)[i]
        end
    end
    println(filters)
    throw(ArgumentError("""No exact match found for $filtername in MIST filterset, please review above list of filters."""))
end
function vegamags(zpt::MISTZeropoints, filter::Union{Symbol, <:AbstractString}, mags)
    nt = zpt(String(filter))
    system = nt.system
    # None of the MIST BCs are in STmags so the only two cases are system == AB or Vega
    if system == "AB"
        return mags .- nt.VegaAB
    else # system == "Vega"
        return mags
    end
end
function abmags(zpt::MISTZeropoints, filter::Union{Symbol, <:AbstractString}, mags)
    nt = zpt(String(filter))
    system = nt.system
    # None of the MIST BCs are in STmags so the only two cases are system == AB or Vega
    if system == "AB"
        return mags
    else # system == "Vega"
        return mags .+ nt.VegaAB
    end
end
function stmags(zpt::MISTZeropoints, filter::Union{Symbol, <:AbstractString}, mags)
    nt = zpt(String(filter))
    system = nt.system
    # None of the MIST BCs are in STmags
    if system == "AB"
        # First convert to Vega, then to ST
        return mags .- nt.VegaAB .+ nt.VegaST
    else # system == "Vega"
        return mags .+ nt.VegaST
    end
end
Mbol(::MISTZeropoints) = 4.74
Lbol(::MISTZeropoints) = 3.828e33

"""
This constant is an instance of [`MISTZeropoints`](@ref BolometricCorrections.MIST.MISTZeropoints).
See the docs for more informations on supported operations.
This constant is returned when calling [`zeropoints`](@ref) on instances of [`MIST.MISTBCGrid`](@ref)
and [`MIST.MISTBCTable`](@ref).

Each set of bolometric corrections is specified on either the Vega or AB magnitude system.
For ease of conversion amongst the AB, ST, and Vega systems this table is provided containing
one line for each filter in our collection. The contents of the file can be used to convert
amongst AB, ST, and Vega magnitude systems as follows.

| filter  | system | mag(Vega/ST) | mag(Vega/AB) |
|---------|--------|--------------|--------------|
| WISE_W1 | Vega   | 6.610497     | 2.665543     |

Above is the information for WISE W1, which is tabulated by default in Vega mags
(as noted in column "system"). To convert WISE mags from Vega to AB is a simple operation:
mag(AB) = mag(Vega) + "mag(Vega/AB)" [column 4 in the file].
For example, a star with WISE\\_W1(Vega) = 0.0 would have WISE\\_W1(AB) = 2.66.

To ease compatibility between `Symbol` and `String` representations, the
mag(Vega/ST) and mag(Vega/AB) columns have been simplified to VegaST and VegaAB.

"""
const zpt = MISTZeropoints(CSV.read(joinpath(@__DIR__, "zeropoints.txt"),
                                    Table; # header=1,
                                    header=[:filter, :system, :VegaST, :VegaAB], skipto=2,
                                    delim=' ', ignorerepeated=true))
# const zeropoints = CSV.read(joinpath(@__DIR__, "zeropoints.txt"), Table; header=1, delim=' ', ignorerepeated=true)
@compat public zpt

# Might be nice to have zpt be a namedtuple, something like
# NamedTuple{Tuple(Symbol(i) for i in BolometricCorrections.MIST.zpt.filter)}(NamedTuple{keys(row)[2:end]}(values(row)[2:end]) for row in rows(BolometricCorrections.MIST.zpt))
# """
# All filters available in the MIST BC grid.
# """
# const filters = zpt.filter
# @compat public filters

# Shared utilities: read a post-processed MIST BC file into a Table
read_mist_bc_processed(fname::AbstractString) =
    CSV.read(fname, Table; buffer_in_memory=true)  # decompress in memory
# Given a datadep directory path, return path to the post-processed .gz file
mist_processed_fname(fname::AbstractString) = joinpath(fname, last(splitpath(fname)) * ".gz")

################
# Initialization
################
include("init.jl")

# All v1.2-specific types and constants
include("MIST_v1.jl")
# All v2.5-specific types and constants
include("MIST_v2.jl")

"Unique values for dependent variables in the MIST v1.2 (key `v1`) and v2.5 (key `v2`) bolometric correction grids."
const gridinfo = (v1 = gridinfov1, v2 = gridinfov2)
@compat public gridinfo

#################################
# Deprecated aliases for v1 types
#################################

@deprecate MISTChemistry() MISTChemistryv1()
@deprecate MISTBCGrid(grid::AbstractString) MISTBCGridv1(grid)
@deprecate MISTBCTable(grid::MISTBCGridv1, feh, Av) MISTBCTablev1(grid, feh, Av)



end # submodule
