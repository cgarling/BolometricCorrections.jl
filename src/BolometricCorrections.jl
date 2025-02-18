module BolometricCorrections

using ArgCheck: @argcheck
using Compat: @compat # for @compat public <x>
import CSV
import Tables
import TypedTables: Table, columnnames, columns, getproperties

# exports from top-level module
# const test = "asdf" # This is exported at the module level
# export test

# Interpolation utilities
include("interp.jl")

# Bolometric correction grid API

""" `AbstractBCGrid{T <: Real}` is the abstract supertype for all bolometric correction grids. `T` is the data type to use internally and is returned by `eltype`. """
abstract type AbstractBCGrid{T <: Real} end
Base.eltype(::AbstractBCGrid{T}) where T = T
"""
    extrema(grid::AbstractBCGrid)

Returns a `NamedTuple` containing the bounds of the dependent variables in the bolometric correction grids (e.g., [Fe/H], Av).
"""
function Base.extrema(::AbstractBCGrid) end
"""
    Table(grid::AbstractBCGrid)

Returns a `TypedTables.Table` containing the data underlying the bolometric correction grid.
"""
function Table(::AbstractBCGrid) end
"""
    columnnames(grid::AbstractBCGrid)

Returns the column names in the provided bolometric correction grid. Includes dependent variables. See `filternames` to retrieve only the names of the contained photometric filters.
    """
columnnames(grid::AbstractBCGrid) = columnnames(Table(grid))
"""
    columns(grid::AbstractBCGrid)

Returns the columns of the table underlying the bolometric correction grid.
"""
columns(grid::AbstractBCGrid) = columns(Table(grid))
getproperties(grid::AbstractBCGrid, names::Tuple{Vararg{Symbol}}) = getproperties(Table(grid), names)
"""
    filternames(grid::AbstractBCGrid)

Returns a `Vector{String}` containing the names of the photometric filters contained in the provided bolometric correction grid. See `columnnames` if you also want to retrieve names of dependent variable columns.
"""
function filternames(::AbstractBCGrid) end

#################################
# Bolometric correction table API

""" `AbstractBCTable{T <: Real}` is the abstract supertype for all bolometric correction tables with extraneous dependent variables (e.g., [Fe/H], Av) fixed -- should typically only have  dependent variables `logg` and `Teff` remaining. `T` is the data type to use internally and is returned by `eltype`. """
abstract type AbstractBCTable{T <: Real} end
Base.eltype(::AbstractBCTable{T}) where T = T
"""
    (table::AbstractBCTable{T})([::Type{TypedTables.Table},]
                                Teff::AbstractArray{S},
                                logg::AbstractArray{V}) where {T, S <: Real, V <: Real}

All concrete subtypes of `AbstractBCTable` must be callable with `(Teff, logg)` arguments and return the interpolated BCs at those values. This method generalizes that to vectors for `Teff` and `logg` and formats the result into a stacked matrix or a `TypedTables.Table`. The two-argument version that returns a matrix is considerably more efficient.
"""
function (table::AbstractBCTable{T})(Teff::AbstractArray{S}, logg::AbstractArray{V}) where {T, S <: Real, V <: Real}
    @argcheck axes(Teff) == axes(logg)
    @argcheck length(Teff) != 0
    filters = filternames(table)
    U = promote_type(T, S, V)
    if isbitstype(U)
        return reshape(reinterpret(U, [table(Teff[i], logg[i]) for i in eachindex(Teff, logg)]), length(Teff), length(filters))
        # return reshape(reinterpret(promote_type(T, S, V), [table(Teff[i], logg[i]) for i in eachindex(Teff, logg)]), length(filters), length(Teff))
    else
        return reduce(hcat, [table(Teff[i], logg[i]) for i in eachindex(Teff, logg)])
        # This hcat's into SMatrix which we don't want
        # reduce(hcat, table(Teff[i], logg[i]) for i in eachindex(Teff, logg)) 
    end
end
function (table::AbstractBCTable{T})(::Type{Table},
                                     Teff::AbstractArray{S},
                                     logg::AbstractArray{V}) where {T, S <: Real, V <: Real}
    filters = filternames(table)
    return Table(Tables.table(table(Teff, logg); header=filters))
    # This is correct but the namedtuple creation is inefficient
    # return Table(NamedTuple{filters}(table(Teff[i], logg[i])) for i in eachindex(Teff, logg))
    # return Table(NamedTuple{filters,NTuple{length(filters),T}}(table(Teff[i], logg[i])) for i in eachindex(Teff, logg))
end
# function (table::AbstractBCTable{T})(::Type{DataFrame},
#                                      Teff::AbstractArray{S},
#                                      logg::AbstractArray{V}) where {T, S <: Real, V <: Real}
#     filters = filternames(table)
#     return DataFrame(table(Teff, logg), SVector(filters))
# end
"""
    extrema(table::AbstractBCTable)

Returns a `NamedTuple` containing the bounds of the dependent variables in the bolometric correction table (e.g., `logg`, `Teff`).
"""
function Base.extrema(::AbstractBCTable) end
"""
    Table(table::AbstractBCTable)

Returns a `TypedTables.Table` containing the data underlying the bolometric correction table.
"""
function Table(::AbstractBCTable) end
"""
    columnnames(table::AbstractBCTable)

Returns the column names in the provided bolometric correction table.
"""
columnnames(table::AbstractBCTable) = columnnames(Table(table))
"""
    columns(table::AbstractBCTable)

Returns the columns of the table underlying the bolometric correction table.
"""
columns(table::AbstractBCTable) = columns(Table(table))
getproperties(table::AbstractBCTable, names::Tuple{Vararg{Symbol}}) = getproperties(Table(table), names)
"""
    filternames(table::AbstractBCGrid)

Returns a `NTuple{N, Symbol}` containing the names of the photometric filters contained in the provided bolometric correction table. See `columnnames` if you also want to retrieve names of dependent variable columns.
"""
function filternames(::AbstractBCTable) end

# Top-level API exports
export Table, columnnames, columns, getproperties, filternames

# Include submodules
include("YBC/YBC.jl")
using .YBC
# exports from YBC
include("MIST/MIST.jl")
using .MIST
@compat public MIST
export MISTBCGrid, MISTBCTable
# exports from MIST



end # module
