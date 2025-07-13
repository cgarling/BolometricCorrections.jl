"""
YBC submodule exposing the bolometric corrections computed by the YBC team based on the white dwarf libraries of [Koester2010](@cite) and [Tremblay2009](@cite) which used Gaia DR2 parallaxes to obtain absolute photometry. As the properties are not expected to vary as a function of initial stellar metallicity, a single table is provided.
"""
module KoesterWD

using ...BolometricCorrections: repack_submatrix, AbstractBCTable, AbstractBCGrid, interp1d
import ...BolometricCorrections: zeropoints, filternames, chemistry # , Y_p, X, X_phot, Y, Y_phot, Z, Z_phot, MH # vegamags, abmags, stmags, Mbol, Lbol
using ..YBC: dtype, pull_table, parse_filterinfo, check_prefix

using ArgCheck: @argcheck
using Compat: @compat
using FITSIO: FITS
using Interpolations: cubic_spline_interpolation, Throw, Flat
using StaticArrays: SVector

# export ...

""" `NTuple{4, Symbol}` listing the dependent variables in the YBC.KoesterWD BC grid. """
const _dependents = (:logTeff, :logg, :Av, :Rv)
""" A_v values in files. For each filter "J", each fits file will have columns "J", "J_Av0.5", "J_Av1", and so on."""
const _Av = (0, 0.5, 1, 2, 5, 10, 20) # Mix of float and integer makes parsing FITS columns easier later
const _logTeff = range(convert(dtype, 3.7), convert(dtype, 4.9); step=convert(dtype, 0.01))
const _logg = range(convert(dtype, 6.5), convert(dtype, 9.5); step=convert(dtype, 0.25))

"""Unique values for dependent variables in the YBC.KoesterWD bolometric correction grid."""
const gridinfo = (logTeff = _logTeff,
                  logg = _logg,
                  Av = _Av,
                  Rv = dtype[3.1])
@compat public gridinfo

"""
    check_vals(Av, gridinfo::NamedTuple)

Validate that ``A_V`` value `Av` is valid for the YBC.KoesterWD white dwarf BC table.
This differs from the main YBC.check_vals implementation in that YBC.KoesterWD does not
vary as a function of metallicity and so we do not check the metallicity value.

```jldoctest
julia> using BolometricCorrections.YBC.KoesterWD: check_vals, gridinfo

julia> check_vals(0.0, gridinfo) # Check passes, returns nothing

julia> using Test: @test_throws, Pass

julia> @test_throws(ArgumentError, check_vals(100.0, gridinfo)) isa Pass # Invalid `Av`, throws error
true
```
"""
function check_vals(Av, gridinfo::NamedTuple)
    Av_ext = extrema(gridinfo.Av)
    if Av < first(Av_ext) || Av > last(Av_ext)
        throw(ArgumentError("Provided A_v $Av is outside the bounds for the BC grid $Av_ext"))
    end
end

#########################################################

"""
    KoesterWDYBCGrid(grid::AbstractString) <: AbstractBCGrid

Load and return the YBC Koester white dwarf bolometric corrections for the given photometric system `grid`,
which must be a valid entry in `BolometricCorrections.YBC.systems`.
This type is used to create instances of 
[`KoesterWDYBCTable`](@ref BolometricCorrections.YBC.KoesterWD.KoesterWDYBCTable)
that have fixed dependent grid variables (Av). This can be done either by calling an instance of
`KoesterWDYBCGrid` with an `Av` argument or by using the appropriate constructor for
[`KoesterWDYBCTable`](@ref BolometricCorrections.YBC.KoesterWD.KoesterWDYBCTable).
Note that there is no variation with metallicity
for this grid as it is specialized for white dwarfs only.

```jldoctest
julia> using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid

julia> grid = KoesterWDYBCGrid("acs_wfc")
YBC Koester white dwarf bolometric correction grid for photometric system YBC/acs_wfc.

julia> grid(0.11) # Can be called to construct table with interpolated Av
YBC Koester white dwarf bolometric correction table with for system YBC/acs_wfc with V-band extinction 0.11
```
"""
struct KoesterWDYBCGrid{A <: Number, C <: AbstractVector{A}, N} <: AbstractBCGrid{A}
    data::Vector{Matrix{A}} # A should be Float32
    mag_zpt::C
    systems::Vector{String}
    name::String
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end

function KoesterWDYBCGrid(data::Vector{Matrix{A}}, mag_zpt::AbstractArray{<:Number}, systems, name, filternames) where {A}
    return KoesterWDYBCGrid(data, convert.(A, mag_zpt), String.(systems), String(name), tuple(Symbol.(filternames)...))
end

function KoesterWDYBCGrid(grid::AbstractString; prefix::AbstractString="YBC")
    check_prefix(prefix)
    path = pull_table(String(grid), String(prefix))
    files = filter(x->occursin("Koester", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    elseif length(files) != 1
        error("Number of files found for grid $grid is $(length(files)), expected 1. Data may be corrupted. \\
        Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Koester WD only has one file
    file = first(files)
    # Read all data and pack into dense matrix
    data = Vector{Matrix{dtype}}(undef, length(gridinfo.Av))
    FITS(file, "r") do f
        for i in eachindex(gridinfo.Av)
            Av = gridinfo.Av[i]
            # Figure out FITS file column name for given Av
            Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(Av)
            data[i] = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
        end
    end
    return KoesterWDYBCGrid(data, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, filternames)
end
(grid::KoesterWDYBCGrid)(Av::Real) = KoesterWDYBCTable(grid, Av)
Base.show(io::IO, z::KoesterWDYBCGrid) = print(io, "YBC Koester white dwarf bolometric correction grid for photometric system $(z.name).")
# function Table(grid::KoesterWDYBCGrid)
#     data = grid.data
#     tables = Vector{Table}(undef, length(data))
# end
Base.extrema(::KoesterWDYBCGrid) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                    logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                    Av = (first(gridinfo.Av), last(gridinfo.Av)),
                                    Rv = (first(gridinfo.Rv), last(gridinfo.Rv)))
filternames(grid::KoesterWDYBCGrid) = grid.filters
chemistry(::Type{<:KoesterWDYBCGrid}) = missing
# zeropoints(::KoesterWDYBCGrid) = zpt

#########################################################
# A single BC table, with fixed [M/H] and Av

"""
    KoesterWDYBCTable(grid::AbstractString, Av::Real) <: AbstractBCTable

Interpolates the YBC bolometric corrections based on the white dwarf
libraries of [Koester2010](@citet) and [Tremblay2009](@citet) which used 
Gaia DR2 parallaxes to obtain absolute photometry to a fixed V-band 
extinction (`Av`), leaving only `Teff` and `logg` as dependent
variables (only one value of `Rv` is provided). Returns an instance 
that is callable with arguments `(Teff [K], logg [cgs])` to interpolate 
the bolometric corrections as a function
of temperature and surface gravity.

```jldoctest
julia> using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid, KoesterWDYBCTable

julia> grid = KoesterWDYBCGrid("acs_wfc")
YBC Koester white dwarf bolometric correction grid for photometric system YBC/acs_wfc.

julia> table = KoesterWDYBCTable(grid, 0.011) # Interpolate table from full grid
YBC Koester white dwarf bolometric correction table with for system YBC/acs_wfc with V-band extinction 0.011

julia> length(table(10_250, 7.65)) == 12 # Returns BC in each filter
true

julia> size(table([10_250, 10_450], [7.65, 7.75])) # `table(array, array)` is also supported
(12, 2)

julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

julia> table(Table, [10_250, 10_450], [7.65, 7.75]) isa Table
true

julia> KoesterWDYBCTable("acs_wfc", 0.5) isa KoesterWDYBCTable # Can construct without a KoesterWDYBCGrid
true
```
"""
struct KoesterWDYBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    Av::A
    mag_zpt::Vector{A}
    systems::Vector{String}
    name::String
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
function KoesterWDYBCTable(Av::Real, mag_zpt::Vector{<:Real}, systems, name, itp, filters)
    T = dtype # promote_type(typeof(Av), eltype(mag_zpt))
    return KoesterWDYBCTable(convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), itp, filters)
end
Base.show(io::IO, z::KoesterWDYBCTable) = print(io, "YBC Koester white dwarf bolometric correction table with for system $(z.name) with V-band extinction ", z.Av)
filternames(table::KoesterWDYBCTable) = table.filters
chemistry(::Type{<:KoesterWDYBCTable}) = missing
# zeropoints(table::KoesterWDYBCTable) = table.mag_zpt

# Interpolations uses `bounds` to return interpolation domain
# We will just use the hard-coded grid bounds; extremely fast
Base.extrema(::KoesterWDYBCTable) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                     logg = (first(gridinfo.logg), last(gridinfo.logg)))
(table::KoesterWDYBCTable)(Teff::Real, logg::Real) = table.itp(logg, log10(Teff))
# Data are naturally Float32 -- convert Float64 args for faster evaluation (~35% faster)
(table::KoesterWDYBCTable)(Teff::Float64, logg::Float64) = table(convert(dtype, Teff), convert(dtype, logg))
# to broadcast over both teff and logg, you do table.(teff, logg')

function KoesterWDYBCTable(grid::AbstractString, Av::Real; prefix::AbstractString="YBC")
    grid, prefix = String(grid), String(prefix)
    check_prefix(prefix)
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use KoesterWDYBCGrid for grid interpolation."
    path = pull_table(grid, prefix)
    files = filter(x->occursin("Koester", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    elseif length(files) != 1
        error("Number of files found for grid $grid is $(length(files)), expected 1. Data may be corrupted. \\
        Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Koester WD only has one file
    goodfile = first(files)
    # Figure out FITS file column name for given Av
    Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
    # Access FITS file
    FITS(goodfile, "r") do f
        data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
        # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
        newdata = repack_submatrix(data, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filternames)))
        itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Throw())
        return KoesterWDYBCTable(Av, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
    end
end

function KoesterWDYBCTable(grid::KoesterWDYBCGrid, Av::Real)
    check_vals(Av, gridinfo)
    filters = filternames(grid)
    data = grid.data
    Av_vec = SVector(gridinfo.Av) # Need vector to use searchsortedfirst
    if Av ∈ gridinfo.Av
        # Exact values are in grid; no interpolation necessary
        submatrix = data[searchsortedfirst(Av_vec, Av)]
    else
        Av_idx = searchsortedfirst(Av_vec, Av) - 1
        mat1 = data[Av_idx]
        mat2 = data[Av_idx + 1]
        submatrix = interp1d(Av, Av_vec[Av_idx], Av_vec[Av_idx + 1], mat1, mat2)       
    end
    newdata = repack_submatrix(submatrix, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filters)))
    itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
    return KoesterWDYBCTable(Av, grid.mag_zpt, grid.systems, grid.name, itp, filters)
end
KoesterWDYBCTable(grid::KoesterWDYBCGrid, Av::Float64) = KoesterWDYBCTable(grid, convert(dtype, Av))

end # module