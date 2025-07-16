"""
YBC submodule exposing the bolometric corrections computed by the YBC team based on the ATLAS9 atmospheres [Kurucz2014](@citet) presented in [Castelli2003](@citet). These models assume [Grevesse1998](@citet) solar chemical abundances.
"""
module ATLAS9

using ArgCheck: @argcheck
using Compat: @compat
using FITSIO: FITS
using Interpolations: cubic_spline_interpolation, Throw, Flat
using StaticArrays: SVector


using ...BolometricCorrections: repack_submatrix, AbstractBCTable, AbstractBCGrid, interp1d, interp2d, AbstractChemicalMixture
import ...BolometricCorrections: zeropoints, filternames, chemistry, Y_p, X, X_phot, Y, Y_phot, Z, Z_phot, MH # vegamags, abmags, stmags, Mbol, Lbol
using ..YBC: HardwareNumeric, dtype, pull_table, parse_filterinfo, check_prefix, check_vals

export ATLAS9YBCTable, ATLAS9YBCGrid, ATLAS9Chemistry

""" `NTuple{5, Symbol}` listing the dependent variables in the YBC.ATLAS9 BC grid. """
const _dependents = (:logTeff, :logg, :MH, :Av, :Rv)
""" A_v values in files. For each filter "J", each fits file will have columns "J", "J_Av0.5", "J_Av1", and so on."""
const _Av = (0, 0.5, 1, 2, 5, 10, 20) # Mix of float and integer makes parsing FITS columns easier later
""" Unique values of [M/H] for the ATLAS9 models. """
const _mh = dtype[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5]
const _logTeff = range(convert(dtype, 3.54), convert(dtype, 4.7); step=convert(dtype, 0.01))
const _logg = range(convert(dtype, 0), convert(dtype, 5.0); step=convert(dtype, 0.5))

"""Unique values for dependent variables in the YBC.ATLAS9 bolometric correction grid."""
const gridinfo = (logTeff = _logTeff,
                  logg = _logg,
                  MH = _mh,
                  Av = _Av,
                  Rv = dtype[3.1])
@compat public gridinfo

"""
    _parse_filename(f::AbstractString)

Return [M/H] of YBC ATLAS9 model given a filename (example: "Avodonnell94Rv3.1fm05k2odfnew.BC.fits").

Note that the YBC repository says
`fm05k2odfnew`: `m05` means [M/H]=-0.05 (`p05` means [M/H]=0.05).

But that would make the dynamic range of [M/H] very low, and if you go to the [ATLAS9 page](https://wwwuser.oats.inaf.it/fiorella.castelli/grids.html)
linked on the repo, they list file formats like `am20k2c125odfnew [-2.0], vturb=2.0 km/s models` indicating
`m20` should be `[M/H] = -2.0` which makes more sense.
"""
function _parse_filename(f::AbstractString)
    f = basename(f)
    @argcheck occursin("odfnew", f)
    mh_str = split(split(f, "3.1f")[2], "k2")[1]
    sign = mh_str[1] == 'p' ? 1 : (mh_str[1] == 'm' ? -1 : error("malformed mh_str $mh_str"))
    mh = parse(dtype, mh_str[2:end]) / 10 * sign
    return (MH = mh,)
end

#########################################################

"""
    ATLAS9YBCGrid(grid::AbstractString)

Load and return the YBC ATLAS9 bolometric corrections for the given photometric system `grid`,
which must be a valid entry in `BolometricCorrections.YBC.systems`.
This type is used to create instances of [`ATLAS9YBCTable`](@ref) that have fixed dependent
grid variables (\\[M/H\\], Av). This can be done either by calling an instance of
`ATLAS9YBCGrid` with `(mh, Av)` arguments or by using the appropriate constructor for [`ATLAS9YBCTable`](@ref).

```jldoctest
julia> grid = ATLAS9YBCGrid("acs_wfc")
YBC ATLAS9 bolometric correction grid for photometric system YBC/acs_wfc.

julia> grid(-1.01, 0.11) # Can be called to construct table with interpolated [M/H], Av
YBC ATLAS9 bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.11
```
"""
struct ATLAS9YBCGrid{A <: Number, C <: AbstractVector{A}, N} <: AbstractBCGrid{A}
    data::Matrix{Matrix{A}} # A should be Float32
    mag_zpt::C
    systems::Vector{String}
    name::String
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end

function ATLAS9YBCGrid(data::Matrix{Matrix{A}}, mag_zpt::AbstractArray{<:Number}, systems, name, filternames) where {A}
    return ATLAS9YBCGrid(data, convert.(A, mag_zpt), String.(systems), String(name), tuple(Symbol.(filternames)...))
end

function ATLAS9YBCGrid(grid::AbstractString; prefix::AbstractString="YBC")
    check_prefix(prefix)
    path = pull_table(String(grid), String(prefix))
    files = filter(x->occursin("odfnew", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    elseif length(files) != length(gridinfo.MH)
        error("Number of files found for grid $grid is $(length(files)), expected $(length(gridinfo.MH)). Data may be corrupted. \\
        Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Sort files by [M/H] value
    idxs = sortperm([_parse_filename(file).MH for file in files])
    files = files[idxs]
    # Check that MH values as parsed from filenames match gridinfo.MH
    if ~mapreduce(isapprox, &, gridinfo.MH, _parse_filename(file).MH for file in files)
        error("File [M/H] values not as expected -- please report.")
    end

    # Read all data and pack into dense matrix
    data = Matrix{Matrix{dtype}}(undef, length(files), length(gridinfo.Av))
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
    return ATLAS9YBCGrid(data, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, filternames)
end
(grid::ATLAS9YBCGrid)(mh::Real, Av::Real) = ATLAS9YBCTable(grid, mh, Av)
Base.show(io::IO, z::ATLAS9YBCGrid) = print(io, "YBC ATLAS9 bolometric correction grid for photometric system $(z.name).")
# function Table(grid::ATLAS9YBCGrid)
#     data = grid.data
#     tables = Vector{Table}(undef, length(data))
# end
Base.extrema(::ATLAS9YBCGrid) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                  logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                  MH = (first(gridinfo.MH), last(gridinfo.MH)),
                                  Av = (first(gridinfo.Av), last(gridinfo.Av)),
                                  Rv = (first(gridinfo.Rv), last(gridinfo.Rv)))
filternames(grid::ATLAS9YBCGrid) = grid.filters
# zeropoints(::ATLAS9YBCGrid) = zpt

#########################################################
# A single BC table, with fixed [M/H] and Av

"""
    ATLAS9YBCTable(grid::ATLAS9YBCGrid, mh::Real, Av::Real)

Interpolates the YBC ATLAS9 bolometric corrections in `grid` to a fixed value of \\[M/H\\]
(`mh`) and V-band extinction (`Av`), leaving only `Teff` and `logg` as dependent
variables (the YBC ATLAS9 BCs have only one `Rv` value). Returns an instance that is callable
with arguments `(Teff [K], logg [cgs])` to interpolate the bolometric corrections as a function
of temperature and surface gravity.

    ATLAS9YBCTable(grid::AbstractString, mh::Real, Av::Real)

Loads the data necessary to construct the BC table for the provided `grid` (e.g., `"acs_wfc"`) 
at \\[M/H\\] = `mh` and V-band extinction `Av`. This method does not support interpolation 
in metallicity or extinction, so the arguments `mh` and `Av` must be among the values 
provided by ATLAS9 (see `BolometricCorrections.YBC.ATLAS9.gridinfo.MH`).

```jldoctest
julia> grid = ATLAS9YBCGrid("acs_wfc")
YBC ATLAS9 bolometric correction grid for photometric system YBC/acs_wfc.

julia> table = ATLAS9YBCTable(grid, -1.01, 0.011) # Interpolate table from full grid
YBC ATLAS9 bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.011

julia> length(table(4025, 0.01)) == 12 # Returns BC in each filter
true

julia> size(table([4025, 4225], [0.01, 0.02])) # `table(array, array)` is also supported
(12, 2)

julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

julia> table(Table, [4025, 4225], [0.01, 0.02]) isa Table
true

julia> ATLAS9YBCTable("acs_wfc", -2.0, 0.5) isa ATLAS9YBCTable # Can construct without a ATLAS9YBCGrid
true
```
"""
struct ATLAS9YBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    MH::A
    Av::A
    mag_zpt::Vector{A}
    systems::Vector{String}
    name::String
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
function ATLAS9YBCTable(MH::Real, Av::Real, mag_zpt::Vector{<:Real}, systems, name, itp, filters)
    T = dtype # promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
    return ATLAS9YBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), itp, filters)
end
Base.show(io::IO, z::ATLAS9YBCTable) = print(io, "YBC ATLAS9 bolometric correction table for system $(z.name) with [M/H] ",
                                              z.MH, " and V-band extinction ", z.Av)
filternames(table::ATLAS9YBCTable) = table.filters
# zeropoints(table::ATLAS9YBCTable) = table.mag_zpt

# Interpolations uses `bounds` to return interpolation domain
# We will just use the hard-coded grid bounds; extremely fast
Base.extrema(::ATLAS9YBCTable) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                  logg = (first(gridinfo.logg), last(gridinfo.logg)))
(table::ATLAS9YBCTable)(Teff::Real, logg::Real) = table.itp(logg, log10(Teff))
# Data are naturally Float32 -- convert hardware numeric args for faster evaluation and guarantee Float32 output
(table::ATLAS9YBCTable)(Teff::HardwareNumeric, logg::HardwareNumeric) = table(convert(dtype, Teff), convert(dtype, logg))
# to broadcast over both teff and logg, you do table.(teff, logg')

function ATLAS9YBCTable(grid::AbstractString, mh::Real, Av::Real; prefix::AbstractString="YBC")
    grid, prefix = String(grid), String(prefix)
    check_prefix(prefix)
    @argcheck mapreduce(isapprox(mh), |, gridinfo.MH) "Provided [M/H] $mh not in available values $(gridinfo.MH); use ATLAS9YBCGrid for grid interpolation."
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use ATLAS9YBCGrid for grid interpolation."
    path = pull_table(grid, prefix)
    files = filter(x->occursin("odfnew", x), readdir(joinpath(path, "regrid"); join=true))
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
        data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
        # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
        newdata = repack_submatrix(data, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filternames)))
        itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
        return ATLAS9YBCTable(mh, Av, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
    end
end

function ATLAS9YBCTable(grid::ATLAS9YBCGrid, mh::Real, Av::Real)
    check_vals(mh, Av, gridinfo)
    filters = filternames(grid)
    data = grid.data

    Av_vec = SVector{length(gridinfo.Av), dtype}(gridinfo.Av) # Need vector to use searchsortedfirst
    MH_vec = gridinfo.MH

    if mh ∈ gridinfo.MH && Av ∈ gridinfo.Av
        # Exact values are in grid; no interpolation necessary
        submatrix = data[searchsortedfirst(MH_vec, mh), searchsortedfirst(Av_vec, Av)]
    else
        if mh ∈ gridinfo.MH
            MH_idx = searchsortedfirst(MH_vec, mh)
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            mat1 = data[MH_idx, Av_idx]
            mat2 = data[MH_idx, Av_idx + 1]
            submatrix = interp1d(Av, Av_vec[Av_idx], Av_vec[Av_idx + 1], mat1, mat2)
        elseif Av ∈ gridinfo.Av
            Av_idx = searchsortedfirst(Av_vec, Av)
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mat1 = data[MH_idx, Av_idx]
            mat2 = data[MH_idx + 1, Av_idx]
            submatrix = interp1d(mh, MH_vec[MH_idx], MH_vec[MH_idx + 1], mat1, mat2)
        else
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            Av1, Av2 = Av_vec[Av_idx], Av_vec[Av_idx + 1]
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mh1, mh2 = MH_vec[MH_idx], MH_vec[MH_idx + 1]

            mat1_1 = data[MH_idx, Av_idx]
            mat2_1 = data[MH_idx + 1, Av_idx]
            mat1_2 = data[MH_idx, Av_idx + 1]
            mat2_2 = data[MH_idx + 1, Av_idx + 1]
            # Perform bilinear interpolation
            submatrix = interp2d(mh, Av, mh1, mh2, Av1, Av2, mat1_1, mat2_1, mat1_2, mat2_2)
        end
    end
    newdata = repack_submatrix(submatrix, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filters)))
    itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
    return ATLAS9YBCTable(mh, Av, grid.mag_zpt, grid.systems, grid.name, itp, filters)
end
ATLAS9YBCTable(grid::ATLAS9YBCGrid, mh::HardwareNumeric, Av::HardwareNumeric) = ATLAS9YBCTable(grid, convert(dtype, mh), convert(dtype, Av))


##############################

"""
    ATLAS9Chemistry()
Returns a singleton struct representing the [ATLAS9](@cite Castelli2003) chemical mixture model.
ATLAS9 assumes the [Grevesse1998](@citet) solar abundances for which photospheric abundances are
equal to protostellar abundances -- this is what they mean when they write "the effects of element
migration at the bottom of the convective zone ... are not observed." In ATLAS9, the helium abundance
is not scaled with the metallicity, so only the hydrogen mass fraction and metal mass fraction change
as a function of metal abundance (i.e., dY/dZ = 0). See also [Girardi2007](@citet).

```jldoctest
julia> using BolometricCorrections.YBC.ATLAS9: ATLAS9Chemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = ATLAS9Chemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> X(chem) == X_phot(chem) # photospheric and protostellar abundances equal (assume no diffusion)
true

julia> MH(chem, Z(chem) * 0.1) ≈ -1.00894760736597
true

julia> Z(chem, -1.00894760736597) ≈ Z(chem) * 0.1
true
```
"""
struct ATLAS9Chemistry <: AbstractChemicalMixture end
chemistry(::Type{<:ATLAS9YBCGrid}) = ATLAS9Chemistry()
chemistry(::Type{<:ATLAS9YBCTable}) = ATLAS9Chemistry()
X(::ATLAS9Chemistry) = 0.735
X_phot(::ATLAS9Chemistry) = 0.735
Y(::ATLAS9Chemistry) = 0.248
Y_phot(::ATLAS9Chemistry) = 0.248
Y_p(::ATLAS9Chemistry) = missing # not defined
Z(::ATLAS9Chemistry) = 0.017
Z_phot(::ATLAS9Chemistry) = 0.017
Y(mix::ATLAS9Chemistry, Zval) = Y(mix) # Constant Y with changing Z, dY/dZ = 0
MH(mix::ATLAS9Chemistry, Zval) = log10((Zval / X(mix, Zval)) * (X_phot(mix) / Z_phot(mix)))
Z(mix::ATLAS9Chemistry, MHval) = (1 - Y(mix)) / (X(mix) / (exp10(MHval) * Z(mix)) + 1)

end # module