module BolometricCorrections

import CSV
import TypedTables: Table, columnnames, columns, getproperties

# exports from top-level module
# const test = "asdf" # This is exported at the module level
# export test


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

include("YBC/YBC.jl")
using .YBC
# exports from YBC
include("MIST/MIST.jl")
using .MIST
# exports from MIST



end # module
