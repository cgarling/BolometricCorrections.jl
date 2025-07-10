"""
Submodule enabling interaction with the MIST grid of stellar bolometric corrections.
"""
module MIST

# using ..BolometricCorrections: Table, columnnames # relative path for parent module
using ..BolometricCorrections: AbstractBCGrid, AbstractBCTable, AbstractZeropoints, AbstractChemicalMixture,
                               interp1d, interp2d, repack_submatrix
import ..BolometricCorrections: zeropoints, filternames, vegamags, abmags, stmags, Mbol, Lbol, Y_p, X, X_phot,
                                Y, Y_phot, Z, Z_phot, MH, chemistry

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

export MISTBCGrid, MISTBCTable, MISTChemistry, X, X_phot, Y, Y_phot, Z, Z_phot, Y_p, MH

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
    return CSV.read(fname, Table; buffer_in_memory=true) # decompress in memory
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
    # If grid is in registry, don't enter parsing stage
    if haskey(registry, grid)
        fname = mist_processed_fname(@datadep_str(grid))
    else # Parse `grid` string
        # Normalize string argument (casefold=true makes lowercase)
        grid = normalize(grid; casefold=true)
        find_func = Base.Fix2(occursin, grid)
        if find_func("jwst") # Only contains NIRCam looks like
            fname = mist_processed_fname(datadep"MIST_JWST")
            # Cover multiple HST instruments gracefully
        elseif find_func("hst")
            if find_func("wfpc2")
                fname = mist_processed_fname(datadep"MIST_HST_WFPC2")
            elseif find_func("wfc3")
                fname = mist_processed_fname(datadep"MIST_HST_WFC3")
            elseif mapreduce(find_func, &, ("acs", "hrc"))
                fname = mist_processed_fname(datadep"MIST_HST_ACS_HRC")
            elseif mapreduce(find_func, &, ("acs", "wfc"))
                fname = mist_processed_fname(datadep"MIST_HST_ACS_WFC")
            else # Name not fully specified, error
                throw(ArgumentError("""Requested grid "$grid" unclear.
                                   Supported HST bolometric correction grids are "hst_acs_wfc" for the \
                                   ACS Wide Field Channel, "hst_acs_hrc" for the ACS High-Resolution \
                                   Channel, "hst_wfpc2" for the Wide Field and Planetary Camera 2, \
                                   and "hst_wfc3" for the Wide Field Camera 3."""))
            end
        elseif find_func("lsst")
            fname = mist_processed_fname(datadep"MIST_LSST")
        elseif find_func("wise")
            fname = mist_processed_fname(datadep"MIST_WISE")
        elseif find_func("washington")
            fname = mist_processed_fname(datadep"MIST_Washington")
        elseif find_func("swift")
            fname = mist_processed_fname(datadep"MIST_Swift")
        elseif find_func("panstarrs")
            fname = mist_processed_fname(datadep"MIST_PanSTARRS")
        elseif find_func("megacam")
            fname = mist_processed_fname(datadep"MIST_CFHT_MegaCam")
        elseif find_func("skymapper")
            fname = mist_processed_fname(datadep"MIST_SkyMapper")
        elseif find_func("uvit")
            fname = mist_processed_fname(datadep"MIST_UVIT")
        elseif find_func("ukidss")
            fname = mist_processed_fname(datadep"MIST_UKIDSS")
        elseif find_func("vista")
            fname = mist_processed_fname(datadep"MIST_VISTA")
        elseif mapreduce(find_func, |, ("johnson", "cousins", "bessell", "2mass", "kepler", "hipparcos", "tycho",
                                        "gaia", "tess"))
            fname = mist_processed_fname(datadep"MIST_UBVRIplus")
        elseif find_func("hsc")
            fname = mist_processed_fname(datadep"MIST_HSC")
        elseif find_func("spitzer")
            fname = mist_processed_fname(datadep"MIST_Spitzer")
        elseif mapreduce(find_func, |, ("sdss", "sloan"))
            fname = mist_processed_fname(datadep"MIST_SDSS")
        elseif find_func("galex")
            fname = mist_processed_fname(datadep"MIST_GALEX")
        elseif find_func("iphas")
            fname = mist_processed_fname(datadep"MIST_IPHAS")
        elseif mapreduce(find_func, |, ("splus", "s+"))
            fname = mist_processed_fname(datadep"MIST_SPLUS")
        elseif mapreduce(find_func, |, ("wfirst", "roman")) # based on May 2018 preliminary filter set
            @info """The WFIRST (now the Nancy Grace Roman Telescope) bolometric corrections are based \
                     on a preliminary filter set from 2018."""
            fname = mist_processed_fname(datadep"MIST_WFIRST")
        elseif find_func("decam")
            fname = mist_processed_fname(datadep"MIST_DECam")
        end
    end
    return MISTBCGrid(read_mist_bc_processed(fname), fname)
end
(grid::MISTBCGrid)(feh::Real, Av::Real) = MISTBCTable(grid, feh, Av)
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
zeropoints(::MISTBCGrid) = zpt


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
zeropoints(::MISTBCTable) = zpt
# Interpolations uses `bounds` to return interpolation domain
# Could also just query _mist_Teff and _mist_logg
Base.extrema(table::MISTBCTable) = (Teff = extrema(table.itp.itp.knots[2]), logg = extrema(table.itp.itp.knots[1]))

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

"""
    MISTBCTable(grid::MISTBCGrid, feh::Real, Av::Real)

Interpolates the MIST bolometric corrections in `grid` to a fixed value of \\[Fe/H\\]
(`feh`) and V-band extinction (`Av`), leaving only `Teff` and `logg` as dependent
variables (the MIST BCs have only one `Rv` value). Returns an instance that is callable
with arguments `(Teff, logg)` to interpolate the bolometric corrections as a function
of temperature and surface gravity.

```jldoctest
julia> grid = MISTBCGrid("JWST")
MIST bolometric correction grid for photometric system MIST_JWST

julia> table = MISTBCTable(grid, -1.01, 0.011)
MIST bolometric correction table with [Fe/H] -1.01 and V-band extinction 0.011

julia> length(table(2755, 0.01)) == 29 # Returns BC in each filter
true

julia> size(table([2755, 2756], [0.01, 0.02])) # `table(array, array)` is also supported
(29, 2)

julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

julia> table(Table, [2755, 2756], [0.01, 0.02]) isa Table
true
```
"""
function MISTBCTable(grid::MISTBCGrid, feh::Real, Av::Real)
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
        elseif Av ∈ _mist_Av
            feh_idx = searchsortedfirst(SVector(_mist_feh), feh) - 1
            feh1, feh2 = _mist_feh[feh_idx], _mist_feh[feh_idx + 1]
            # Extract corresponding subtables, retrieve BCs from subtable, convert to matrix
            mat1 = Tables.matrix(getproperties(select_subtable(table, feh1, Av), filters))
            mat2 = Tables.matrix(getproperties(select_subtable(table, feh2, Av), filters))
            # submatrix = @. (mat1 * (feh2 - feh) + mat2 * (feh - feh1)) / (feh2 - feh1)
            submatrix = interp1d(feh, feh1, feh2, mat1, mat2)
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
        end
    end
    # Construct interpolator and return MISTBCTable
    itp = interpolate((SVector(_mist_logg), SVector(_mist_Teff)),
                      repack_submatrix(submatrix, length(_mist_logg), length(_mist_Teff), filters),
                      Gridded(Linear()))
    itp = extrapolate(itp, Flat())
    return MISTBCTable(feh, Av, itp, filters)
end
(table::MISTBCTable)(Teff::Real, logg::Real) = table.itp(logg, Teff)
# to broadcast over both teff and logg, you do table.(teff, logg')

##############################
# Chemical mixture information
"""
    MISTChemistry()
Returns a singleton struct representing the MIST chemical mixture model.
MIST assumes the *protostellar* [Asplund2009](@citet) solar abundances. Sum of protostellar hydrogen, helium,
metal mass fractions from last row of Table 4 sums to 0.9999, not 1 as it should.
To keep calculations consistent, the protostellar values are normalized to sum to 1 here.

```jldoctest
julia> using BolometricCorrections.MIST: MISTChemistry, X, Y, Z, X_phot, Y_phot, Z_phot,
                                         MH;

julia> chem = MISTChemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> MH(chem, Z(chem) * 0.1) ≈ -1.0189881255814277
true

julia> Z(chem, -1.0189881255814277) ≈ Z(chem) * 0.1
true
```
"""
struct MISTChemistry <: AbstractChemicalMixture end
chemistry(::MISTBCGrid) = MISTChemistry()
chemistry(::MISTBCTable) = MISTChemistry()
X(::MISTChemistry) = 0.7154 / 0.9999
X_phot(::MISTChemistry) = 0.7381
Y(::MISTChemistry) = 0.2703 / 0.9999
Y_phot(::MISTChemistry) = 0.2485
Y_p(::MISTChemistry) = 0.249
Z(::MISTChemistry) = 0.0142 / 0.9999
Z_phot(::MISTChemistry) = 0.0134
function Y(mix::MISTChemistry, Zval)
    yp = Y_p(mix)
    return yp + ((Y(mix) - yp) / Z(mix)) * Zval
end
function MH(mix::MISTChemistry, Zval)
    Xval = X(mix, Zval)
    solZ = Z(mix)
    solX = X(mix)
    # return log10(Zval / Xval) - log10(solZ / solX)
    # Fuse expression into one log10 call for efficiency
    return log10((Zval / Xval) * (solX / solZ))
end
function Z(mix::MISTChemistry, MHval)
    # [M/H] = log(Z/X)-log(Z/X)☉ with Z☉ = solz
    # Z/X = exp10( [M/H] + log(Z/X)☉ )
    # X = 1 - Y - Z
    # Y ≈ Y_p + ((Y☉ - Y_p) / Z☉) * Z (see Y above)
    # so X ≈ 1 - (Y_p + ((Y☉ - Y_p) / Z☉) * Z) - Z = 1 - Y_p - (((Y☉ - Y_p) / Z☉) * Z) - Z
    # Substitute into line 2,
    # Z / (1 - Y_p - Z - (((Y☉ - Y_p) / Z☉) * Z)) = exp10( [M/H] + log(Z/X)☉ )
    # Z = (1 - Y_p - Z - (((Y☉ - Y_p) / Z☉) * Z)) * exp10( [M/H] + log(Z/X)☉ )
    # let A = exp10( [M/H] + log(Z/X)☉ )
    # Z = (1 - Y_p - Z - (((Y☉ - Y_p) / Z☉) * Z)) * A
    # Z = A - A * Y_p - A * Z - A * (((Y☉ - Y_p) / Z☉) * Z)
    # Z + A * (((Y☉ - Y_p) / Z☉) * Z) + A * Z = A * (1 - Y_p)
    # Z * (1 + A + A * ((Y☉ - Y_p) / Z☉)) = A * (1 - Y_p)
    # Z = A * (1 - Y_p) / (1 + A + A * ((Y☉ - Y_p) / Z☉))
    solX, solY, solZ, yp = X(mix), Y(mix), Z(mix), Y_p(mix)
    zoverx = exp10(MHval + log10(solZ/solX))
    return zoverx * (1 - yp) / (1 + zoverx + zoverx * ((solY - yp) / solZ))
end

end # submodule
