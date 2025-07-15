module BolometricCorrections

using ArgCheck: @argcheck
using Compat: @compat # for @compat public <x>
using StaticArrays: SVector
import CSV
import Tables
import TypedTables: Table, columnnames, columns, getproperties

# exports from top-level module
# const test = "asdf" # This is exported at the module level
# export test

# Interpolation utilities
include("interp.jl")

"""All basic hardware numeric types supported by Julia."""
const AllHardwareNumeric = Union{Int8, Int16, Int32, Int64, Int128,
    UInt8, UInt16, UInt32, UInt64, UInt128, Float16, Float32, Float64}
without(dtype, union::Union = HardwareNumeric) = Union{filter(t -> t !== dtype, Base.uniontypes(union))...}

#################################
# Bolometric correction grid API

""" `AbstractBCGrid{T <: Real}` is the abstract supertype for all bolometric correction grids. `T` is the data type to use internally and is returned by `eltype`. Generally, concrete subtypes should be callable with population properties (e.g., metallicity, reddening, etc.) to interpolate the full grid to these properties, returning a concrete subtype of [`BolometricCorrections.AbstractBCTable`](@ref). As different grids will have different population properties available (e.g., some support different α-element abundances in addition to total metallicity), the call signature to interpolate the grid is specific for each concrete subtype, which include
 - [`MISTBCGrid`](@ref)
 - [`PHOENIXYBCGrid`](@ref)
 - [`ATLAS9YBCGrid`](@ref)
"""
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
"""
    getproperties(grid::AbstractBCGrid, names::Tuple{Vararg{Symbol}})

Returns properties `names` from the provided bolometric correction grid.
"""
getproperties(grid::AbstractBCGrid, names::Tuple{Vararg{Symbol}}) = getproperties(Table(grid), names)
"""
    filternames(grid::AbstractBCGrid)

Returns a `Vector{String}` containing the names of the photometric filters contained in the provided bolometric correction grid. See `columnnames` if you also want to retrieve names of dependent variable columns.
"""
function filternames(::AbstractBCGrid) end

#################################
# Bolometric correction table API

""" `AbstractBCTable{T <: Real}` is the abstract supertype for all bolometric correction tables with extraneous dependent variables (e.g., [Fe/H], Av) fixed -- should typically only have  dependent variables `logg` and `Teff` remaining. `T` is the data type to use internally and is returned by `eltype`.

    (table::AbstractBCTable{T})([::Type{TypedTables.Table},]
                                Teff::AbstractArray{S},
                                logg::AbstractArray{V}) where {T, S <: Real, V <: Real}

All concrete subtypes of `AbstractBCTable` must be callable with `(Teff, logg)` arguments and return the interpolated BCs at those values. This method broadcasts the operation over arrays of `Teff` and `logg` and formats the result into a stacked matrix or a `TypedTables.Table`. The three-argument version that returns a `Table` has a roughly fixed runtime overhead cost of 3--5 μs to perform the type conversion. 
"""
abstract type AbstractBCTable{T <: Real} end
Base.eltype(::AbstractBCTable{T}) where T = T
function (table::AbstractBCTable{T})(Teff::AbstractArray{S}, logg::AbstractArray{V}) where {T, S <: Real, V <: Real}
    @argcheck axes(Teff) == axes(logg)
    @argcheck length(Teff) != 0
    filters = filternames(table)
    U = promote_type(T, S, V)
    if isbitstype(U)
        return reshape(reinterpret(U, [table(Teff[i], logg[i]) for i in eachindex(Teff, logg)]), length(filters), length(Teff))
    else
        return reduce(hcat, [table(Teff[i], logg[i]) for i in eachindex(Teff, logg)])
        # This hcat's into SMatrix which we don't want
        # reduce(hcat, table(Teff[i], logg[i]) for i in eachindex(Teff, logg)) 
    end
end
function (table::AbstractBCTable)(::Type{Table},
                                  Teff::AbstractArray{<:Real},
                                  logg::AbstractArray{<:Real})
    filters = filternames(table)
    # Mostly fixed ~3--5 μs runtime cost
    return Table(Tables.table(PermutedDimsArray(table(Teff, logg), (2,1)); header=filters))
    # More expensive due to data copy with permutedims
    # return Table(Tables.table(permutedims(table(Teff, logg)); header=filters))
    # This is correct but the cost of namedtuple creation scales poorly with length of argument
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
Base.extrema(::AbstractBCTable)
"""
    Table(table::AbstractBCTable)

Returns a `TypedTables.Table` containing the data underlying the bolometric correction table.
"""
Table(::AbstractBCTable)
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
"""
    getproperties(grid::AbstractBCGrid, names::Tuple{Vararg{Symbol}})

Returns properties `names` from the provided bolometric correction grid.
"""
getproperties(table::AbstractBCTable, names::Tuple{Vararg{Symbol}}) = getproperties(Table(table), names)
"""
    filternames(table::AbstractBCGrid)

Returns a `NTuple{N, Symbol}` containing the names of the photometric filters contained in the provided bolometric correction table. See `columnnames` if you also want to retrieve names of dependent variable columns.
"""
function filternames(::AbstractBCTable) end

#########################################
# Zeropoint definition and conversion API

""" `AbstractZeropoints` is the abstract supertype for information regarding the
photometric zeropoints assumed for a particular grid of bolometric corrections and supports
conversion between systems (AB, Vega, ST). """
abstract type AbstractZeropoints end # {T <: Real} end

"""
    zeropoints(grid::AbstractBCGrid)
    zeropoints(table::AbstractBCTable)

Return the correct concrete instance of [`AbstractZeropoints`](@ref BolometricCorrections.AbstractZeropoints)
for the type of `grid` or `table`.

```jldoctest
julia> zeropoints(MISTBCGrid("JWST")) isa BolometricCorrections.MIST.MISTZeropoints
true
```
"""
function zeropoints(grid::AbstractBCGrid) end
"""
    vegamags(zpt::AbstractZeropoints, filter, mags)

Uses the photometric zeropoint information in `zpt` to convert magnitudes `mags`
in the given `filter` to the Vega magnitude system.
"""
function vegamags(zpt::AbstractZeropoints, filter, mags) end
"""
    abmags(zpt::AbstractZeropoints, filter, mags)

Uses the photometric zeropoint information in `zpt` to convert magnitudes `mags`
in the given `filter` to the AB magnitude system.
"""
function abmags(zpt::AbstractZeropoints, filter, mags) end
"""
    stmags(zpt::AbstractZeropoints, filter, mags)

Uses the photometric zeropoint information in `zpt` to convert magnitudes `mags`
in the given `filter` to the ST magnitude system.
"""
function stmags(zpt::AbstractZeropoints, filter, mags) end
"""
    filternames(zpt::AbstractZeropoints)

Returns a `Vector{String}` containing the names of the photometric filters in the table of zeropoints.
"""
function filternames(zpt::AbstractZeropoints) end
"""
    Mbol(zpt::AbstractZeropoints)

Returns the absolute bolometric magnitude of the Sun assumed in the definition of the BC grid.
"""
function Mbol(zpt::AbstractZeropoints) end
"""
    Lbol(zpt::AbstractZeropoints)

Returns the bolometric luminosity [erg / s] of the Sun assumed in the definition of the BC grid.
"""
function Lbol(zpt::AbstractZeropoints) end

#################################
# Chemical mixture API

""" `AbstractChemicalMixture` is the abstract supertype for information regarding the
chemical mixtures available or assumed for a particular grid of bolometric corrections.
It also supports evaluation and conversion between different chemical conventions
(i.e., conversion between metal mass fraction [`Z`](@ref) and logarithmic metallicity
[\\[M/H\\]](@ref MH)). """
abstract type AbstractChemicalMixture end
Base.Broadcast.broadcastable(t::AbstractChemicalMixture) = Ref(t)

"""
    chemistry(mix::AbstractBCTable)
    chemistry(mix::AbstractBCGrid)
Returns the correct concrete instance of `AbstractChemicalMixture` for the
provided bolometric correction grid or table. This provides a convenient 
programmatic way to obtain this chemical information.

```jldoctest
julia> grid = MISTBCGrid("JWST");

julia> chemistry(grid)
BolometricCorrections.MIST.MISTChemistry()

julia> table = grid(-1.5, 0.03);

julia> chemistry(table)
BolometricCorrections.MIST.MISTChemistry()
```
"""
chemistry(mix::AbstractBCTable) = chemistry(typeof(mix))
chemistry(mix::AbstractBCGrid) = chemistry(typeof(mix))
# ↑ generics that call to chemistry(::Type{<:NewType}) which can often be simple

"""
    X(mix::AbstractChemicalMixture)
Returns the **protostellar** solar hydrogen mass fraction assumed in the provided chemical mixture.
"""
function X(mix::AbstractChemicalMixture) end
"""
    X_phot(mix::AbstractChemicalMixture)
Returns the **photospheric** solar hydrogen mass fraction assumed in the provided chemical mixture.
"""
function X_phot(mix::AbstractChemicalMixture) end
"""
    X(mix::AbstractChemicalMixture, Z)
Returns **protostellar** hydrogen mass fraction given the protostellar metal mass `Z` and the
provided chemical mixture. Generic method returns `1 - Y(mix, Z) - Z`. 
"""
function X(mix::AbstractChemicalMixture, Z)
    Yval = Y(mix, Z)
    return 1 - Yval - Z
end

"""
    Y_p(mix::AbstractChemicalMixture)
Returns the primordial helium mass fraction assumed in the provided chemical mixture.
"""
function Y_p(mix::AbstractChemicalMixture) end
"""
    Y(mix::AbstractChemicalMixture)
Returns the **protostellar** solar helium mass fraction. May use [`Z`] internally.
"""
function Y(mix::AbstractChemicalMixture) end
"""
    Y_phot(mix::AbstractChemicalMixture)
Returns the **photospheric** solar helium mass fraction. May use [`Z_phot`] internally.
"""
function Y_phot(mix::AbstractChemicalMixture) end
"""
    Y(mix::AbstractChemicalMixture, Z)
Returns the **protostellar** helium mass fraction given the protostellar metal mass
fraction `Z`. Only valid for grids in which ``Y`` is a function of ``Z``.
"""
function Y(mix::AbstractChemicalMixture, Z) end

"""
    Z(mix::AbstractChemicalMixture)
Returns the **protostellar** solar metal mass fraction assumed in the provided chemical mixture. 
"""
function Z(mix::AbstractChemicalMixture) end
"""
    Z(mix::AbstractChemicalMixture, MH)
Returns the **protostellar** metal mass fraction given the logarithmic metallicity `MH`
and the provided chemical mixture. 
"""
function Z(mix::AbstractChemicalMixture, MH) end
"""
    Z_phot(mix::AbstractChemicalMixture)
Returns the **photospheric** solar metal mass fraction assumed in the provided chemical mixture.
"""
function Z_phot(mix::AbstractChemicalMixture) end
"""
    MH(mix::AbstractChemicalMixture, Z)
Returns the **protostellar** logarithmic metallicity [M/H] = log10(Z/X) - log10(Z⊙/X⊙)
given the metal mass fraction `Z` and the provided chemical mixture.
"""
function MH(mix::AbstractChemicalMixture, Z) end
# function Z(mix::T) where T <: AbstractChemicalMixture
#     @warn "Requested solar protostellar metal mass for chemical mixture model $T. This model does not have a `Z` method implemented, so we are falling back to the photospheric metal mass fraction `Z_phot(mix)`." maxlog=1
#     return Z_phot(mix) # Fallback for unimplemented Z
# end
# function Z_phot(mix::T) where T <: AbstractChemicalMixture
#     @warn "Requested solar photospheric metal mass for chemical mixture model $T. This model does not have a `Z_phot` method implemented, so we are falling back to the protostellar metal mass fraction `Z(mix)`." maxlog=1
#     return Z(mix) # Fallback for unimplemented Z_phot
# end

#################################
# Top-level API exports
export Table, columnnames, columns, getproperties, filternames, zeropoints, vegamags,
    stmags, abmags, Mbol, Lbol, X, X_phot, Y_p, Y, Y_phot, Z, Z_phot, MH, chemistry

# Include submodules
include(joinpath("MIST", "MIST.jl"))
using .MIST
@compat public MIST
export MISTBCGrid, MISTBCTable

include(joinpath("YBC", "YBC.jl"))
using .YBC
@compat public YBC
export PHOENIXYBCTable, PHOENIXYBCGrid, ATLAS9YBCTable, ATLAS9YBCGrid


end # module
