"""
Submodule enabling interaction with the MIST grid of stellar bolometric corrections.
"""
module MIST

# using ..BolometricCorrections: Table, columnnames # relative path for parent module
using ..BolometricCorrections: AbstractBCGrid, AbstractBCTable, interp1d, interp2d
import ..BolometricCorrections: filternames
using ArgCheck: @argcheck
using CodecXz: XzDecompressorStream # Decompress downloaded BCs
using Compat: @compat # for @compat public <x>
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

export MISTBCGrid, MISTBCTable

""" `NTuple{5, Symbol}` listing the dependent variables in the MIST BC grid. """
const _mist_dependents = (:Teff, :logg, :feh, :Av, :Rv)
""" `NTuple{5, Symbol}` giving the order in which the MIST dependent variables iterate in the post-processed data tables. """
const _mist_dependents_order = (:logg, :Teff, :Av, :Rv, :feh)
""" Unique values of `Teff` in the MIST BC tables. """
const _mist_Teff = (2500.0, 2800.0, 3000.0, 3200.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, 60000.0, 70000.0, 80000.0, 90000.0, 100000.0, 110000.0, 120000.0, 130000.0, 140000.0, 150000.0, 160000.0, 170000.0, 180000.0, 190000.0, 200000.0, 300000.0, 400000.0, 500000.0, 600000.0, 700000.0, 800000.0, 900000.0, 1.0e6) # length=70
""" Unique values of `logg` in the MIST BC tables. """
const _mist_logg = (-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5) # length=26
""" Unique values of `Av` in the MIST BC tables. """
const _mist_Av = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0) # length=13
""" Unique values of `Rv` in the MIST BC tables. """
const _mist_Rv = (3.1,) # length=1
""" Unique values of [Fe/H] in the MIST BC tables. """
const _mist_feh = (-4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75) # length=18

# User-facing information on the grid 
""" Unique values for dependent variables in the MIST bolometric correction grid. """
const gridinfo = (Teff = _mist_Teff,
                  logg = _mist_logg,
                  feh = _mist_feh,
                  Av = _mist_Av,
                  Rv = _mist_Rv)
                  # dependents = _mist_dependents,
                  # dependents_order = _mist_dependents_order,
@compat public gridinfo

"""
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
"""
const zeropoints = CSV.read(joinpath(@__DIR__, "zeropoints.txt"), Table; header=1, delim=' ', ignorerepeated=true)
@compat public zeropoints
# """
# All filters available in the MIST BC grid.
# """
# const filters = zeropoints.filter
# @compat public filters
                  

#############################
# Initialization and datadeps
include("init.jl")

"""
    check_vals(feh, Av)

Validate that [Fe/H] value `feh` and ``A_V`` value `Av` are valid for the MIST BC grid.
Throws `DomainError` if check fails, returns `nothing` if check is successful.
"""
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
"""
    MISTBCGrid(grid::AbstractString)

Load and return the MIST bolometric corrections for the given photometric system `grid`.
This type is used to create instances of [`MISTBCTable`](@ref) that have fixed dependent
grid variables (\\[Fe/H\\], Av, Rv). This can be done either by calling an instance of
`MISTBCGrid` with `(feh, Av)` arguments or by using the appropriate constructor for [`MISTBCTable`](@ref).

Examples
```jldoctest
julia> grid = MISTBCGrid("JWST")
MIST bolometric correction grid for photometric system MIST_JWST

julia> grid(-1.01, 0.11) # Can be called to construct table with interpolated [Fe/H], Av
MIST bolometric correction table with [Fe/H] -1.01 and V-band extinction 0.11
```
"""
struct MISTBCGrid{A,B,C <: AbstractString} <: AbstractBCGrid{A}
    table::B # Usually a TypedTables.Table
    filename::C
    function MISTBCGrid(table::B, filename::C) where {B, C}
        A = Base.promote_eltype(first(table))
        new{A, B, C}(table, filename)
    end
end
function MISTBCGrid(grid::AbstractString)
    if mapreduce(x->occursin(x,grid), |, ("JWST", "jwst"))
        fname = mist_processed_fname(datadep"MIST_JWST")
        return MISTBCGrid(read_mist_bc_processed(fname), fname)
    end
end
(grid::MISTBCGrid)(feh::Real, Av::Real) = MISTBCTable(feh, Av, grid)
Base.show(io::IO, z::MISTBCGrid) = print(io, "MIST bolometric correction grid for photometric system ",
                                         splitpath(splitext(z.filename)[1])[end])
Table(grid::MISTBCGrid) = grid.table
# columnnames(grid::MISTBCGrid) = columnnames(Table(grid))
# columns(grid::MISTBCGrid) = columns(Table(grid))
# getproperties(grid::MISTBCGrid, names::Tuple{Vararg{Symbol}}) = getproperties(Table(grid), names) 
# A function that will extract the dependent variables from a MIST BC grid
extract_dependents(grid::MISTBCGrid) = getproperties(grid, _mist_dependents)
function Base.extrema(grid::MISTBCGrid)
    return NamedTuple{_mist_dependents}(extrema(col) for col in columns(extract_dependents(grid)))
end
# filternames(grid::MISTBCGrid) = [string(name) for name in columnnames(grid)[6:end]]
filternames(grid::MISTBCGrid) = columnnames(grid)[length(_mist_dependents)+1:end]


#########################################################
# A single BC table, with fixed feh and Av
struct MISTBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    feh::A
    Av::A
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
MISTBCTable(feh::Real, Av::Real, itp, filters) = MISTBCTable(promote(feh, Av)..., itp, filters)
Base.show(io::IO, z::MISTBCTable) = print(io, "MIST bolometric correction table with [Fe/H] ",
                                          z.feh, " and V-band extinction ", z.Av)
filternames(table::MISTBCTable) = table.filters
# Interpolations uses `bounds` to return interpolation domain
# Could also just query _mist_Teff and _mist_logg
Base.extrema(table::MISTBCTable) = (Teff = extrema(table.itp.knots[2]), logg = extrema(table.itp.knots[1]))

# Extract a subtable out of table where table.feh == feh and table.Av == Av
function select_subtable(table::Table, feh::Real, Av::Real)
    @argcheck feh ∈ _mist_feh && Av ∈ _mist_Av
    # return filter(row -> (row.feh ≈ feh) && (row.Av ≈ Av), table)
    # We can use the known structure of the data table to prevent having to do a filter at all
    # For each unique Rv (1), Av and feh there will be
    # length(Teff) * length(logg) * length(Rv) rows
    nrows_per_Av = length(_mist_Teff) * length(_mist_logg) * length(_mist_Rv) # 1820
    nrows_per_feh = nrows_per_Av * length(_mist_Av)
    # The order of iteration of the dependent variables is given by _mist_dependents_order
    # We need to identify the correct feh first, then the correct Av within that
    idx1 = 1 + nrows_per_feh * (findfirst(==(feh), _mist_feh) - 1)
    idx2 = idx1 + nrows_per_Av * (findfirst(==(Av), _mist_Av) - 1)
    idx3 = idx2 + nrows_per_Av - 1
    result = table[idx2:idx3]
    # @argcheck all(@. result.feh == feh && result.Av == Av)
    return result
end

# Use statically known size from filters argument to repack submatrix
# into a vector of SVectors to pass into interpolator
function repack_submatrix(submatrix::AbstractArray{T},
                          filters::NTuple{N, Symbol}) where {T, N}
    submatrix = reshape(submatrix,
                        length(_mist_logg),
                        length(_mist_Teff),
                        length(filters))
    return [SVector{N, T}(view(submatrix,i,j,:)) for i=axes(submatrix,1),j=axes(submatrix,2)]
end
"""
    MISTBCTable(feh::Real, Av::Real, grid::MISTBCGrid)

Interpolates the MIST bolometric corrections in `grid` to a fixed value of \\[Fe/H\\]
(`feh`) and V-band extinction (`Av`), leaving only `Teff` and `logg` as dependent
variables (the MIST BCs have only one `Rv` value). Returns an instance that is callable
with arguments `(Teff, logg)` to interpolate the bolometric corrections as a function
of temperature and surface gravity.

Examples
```jldoctest
julia> grid = MISTBCGrid("JWST")
MIST bolometric correction grid for photometric system MIST_JWST

julia> table = MISTBCTable(-1.01, 0.011, grid)
MIST bolometric correction table with [Fe/H] -1.01 and V-band extinction 0.011

julia> length(table(2755, 0.01)) == 29 # Returns BC in each filter
true
```
"""
function MISTBCTable(feh::Real, Av::Real, grid::MISTBCGrid)
    check_vals(feh, Av)
    filters = filternames(grid)
    table = Table(grid)

    # Exact values are in grid; no interpolation necessary
    if feh ∈ _mist_feh && Av ∈ _mist_Av
        # Basically all the time is spent in this filter ... 
        # subtable = filter(row -> (row.feh ≈ feh) && (row.Av ≈ Av), table)
        subtable = select_subtable(table, feh, Av)
        # Need to repack the subtable into the correct format for interpolation
        submatrix = Tables.matrix(getproperties(subtable, filters))
        itp = interpolate((SVector(_mist_logg), SVector(_mist_Teff)),
                          repack_submatrix(submatrix, filters),
                          Gridded(Linear()))
        return MISTBCTable(feh, Av, itp, filters)
    else
        # Need to interpolate table to correct [Fe/H] and Av
        # Doing basic bilinear interpolation on the matrix representation
        # of the table to get a new matrix that we can reshape back into a subtable
        if feh ∈ _mist_feh
            Av_idx = searchsortedfirst(SVector(_mist_Av), Av) - 1
            Av1, Av2 = _mist_Av[Av_idx], _mist_Av[Av_idx + 1]
            # Extract corresponding subtables, retrieve BCs from subtable, convert to matrix
            mat1 = Tables.matrix(getproperties(select_subtable(table, feh, Av1), filters))
            mat2 = Tables.matrix(getproperties(select_subtable(table, feh, Av2), filters))
            # submatrix = @. mat1 * (Av2 - Av) / (Av2 - Av1) + mat2 * (Av - Av1) / (Av2 - Av1) # Equivalent, bit slower
            # submatrix = @. (mat1 * (Av2 - Av) + mat2 * (Av - Av1)) / (Av2 - Av1)
            submatrix = interp1d(Av, Av1, Av2, mat1, mat2)
            itp = interpolate((SVector(_mist_logg), SVector(_mist_Teff)),
                              repack_submatrix(submatrix, filters),
                              Gridded(Linear()))
            return MISTBCTable(feh, Av, itp, filters)
        elseif Av ∈ _mist_Av
            feh_idx = searchsortedfirst(SVector(_mist_feh), feh) - 1
            feh1, feh2 = _mist_feh[feh_idx], _mist_feh[feh_idx + 1]
            # Extract corresponding subtables, retrieve BCs from subtable, convert to matrix
            mat1 = Tables.matrix(getproperties(select_subtable(table, feh1, Av), filters))
            mat2 = Tables.matrix(getproperties(select_subtable(table, feh2, Av), filters))
            # submatrix = @. (mat1 * (feh2 - feh) + mat2 * (feh - feh1)) / (feh2 - feh1)
            submatrix = interp1d(feh, feh1, feh2, mat1, mat2)
            itp = interpolate((SVector(_mist_logg), SVector(_mist_Teff)),
                              repack_submatrix(submatrix, filters),
                              Gridded(Linear()))
            return MISTBCTable(feh, Av, itp, filters)
        else
            # Identify nearest [Fe/H] and Av values in the table;
            # all dependent variables are interated in sorted order
            # in the post-processed data tables
            feh_idx = searchsortedfirst(SVector(_mist_feh), feh) - 1
            feh1, feh2 = _mist_feh[feh_idx], _mist_feh[feh_idx + 1]
            Av_idx = searchsortedfirst(SVector(_mist_Av), Av) - 1
            Av1, Av2 = _mist_Av[Av_idx], _mist_Av[Av_idx + 1]
            # Extract corresponding subtables, retrieve BCs from subtable, convert to matrix
            mat1_1 = Tables.matrix(getproperties(select_subtable(table, feh1, Av1), filters))
            mat2_1 = Tables.matrix(getproperties(select_subtable(table, feh2, Av1), filters))
            mat1_2 = Tables.matrix(getproperties(select_subtable(table, feh1, Av2), filters))
            mat2_2 = Tables.matrix(getproperties(select_subtable(table, feh2, Av2), filters))
            # Perform bilinear interpolation
            submatrix = interp2d(feh, Av, feh1, feh2, Av1, Av2, mat1_1, mat2_1, mat1_2, mat2_2)

            # # This approach ~550 μs, 350μs for interp2d
            # # Interpolations.jl expects for grid points (x,y), z has shape (length(x), length(y))
            # bigmat = Matrix{Matrix{Float64}}(undef, 2, 2)
            # bigmat[1,1] = mat1_1
            # bigmat[2,1] = mat2_1
            # bigmat[1,2] = mat1_2
            # bigmat[2,2] = mat2_2
            # submatrix = interpolate((SVector(feh1, feh2), SVector(Av1, Av2)),
            #                         # [mat1_1, mat1_2, mat2_1, mat2_2],
            #                         bigmat,
            #                         Gridded(Linear()))(feh, Av)

            # Construct interpolator and return MISTBCTable
            itp = interpolate((SVector(_mist_logg), SVector(_mist_Teff)),
                              repack_submatrix(submatrix, filters),
                              Gridded(Linear()))
            return MISTBCTable(feh, Av, itp, filters)
        end
    end
end
(table::MISTBCTable)(Teff::Real, logg::Real) = table.itp(logg, Teff)
# to broadcast over both teff and logg, you do table.(teff, logg')

end # submodule
