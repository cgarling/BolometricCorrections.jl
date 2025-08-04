
"""
YBC submodule exposing the bolometric corrections based on PHOENIX BT-Settl atmospheres computed by the YBC team. These models assume [Asplund2009](@citet) solar chemical abundances.

The main reference article for these models is [Allard2012](@citet). [Allard2013](@citet) discusses BT-Settl models with the solar composition from [Caffau2011](@citet), but this module uses the model with the [Asplund2009](@citet) abundances.
"""
module PHOENIX

using ...BolometricCorrections: repack_submatrix, AbstractBCTable, AbstractBCGrid, interp1d, interp2d
import ...BolometricCorrections: zeropoints, filternames, chemistry, Z, MH, gridname # Y_p, X, X_phot, Y, Y_phot, Z_phot, vegamags, abmags, stmags, Mbol, Lbol
using ...BolometricCorrections.MIST: MISTChemistry # MIST and YBC PHOENIX both use Asplund2009 abundances, so just use MISTChemistry
using ..YBC: HardwareNumeric, dtype, pull_table, parse_filterinfo, check_prefix, check_vals

using ArgCheck: @argcheck
using Compat: @compat
using FITSIO: FITS, read_header, colnames
# using Interpolations: interpolate, extrapolate, Flat, Throw, BSpline, Cubic, Line, OnGrid
using Interpolations: cubic_spline_interpolation, Throw, Flat
# import CSV
using Printf: @sprintf # Formatted conversion of floats to strings
using StaticArrays: SVector
using TypedTables: Table

export PHOENIXYBCTable, PHOENIXYBCGrid

""" `NTuple{5, Symbol}` listing the dependent variables in the YBC.PHOENIX BC grid. """
const _dependents = (:logTeff, :logg, :MH, :Av, :Rv)
""" A_v values in files. For each filter "J", each fits file will have columns "J", "J_Av0.5", "J_Av1", and so on."""
const _Av = (0, 0.5, 1, 2, 5, 10, 20) # Mix of float and integer makes parsing FITS columns easier later
# const _Av = ["0.5", "1", "2", "5", "10", "20"]
""" Unique values of [M/H] for the PHOENIX BT-Settl models. """
const _mh = dtype[-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.0, 0.3, 0.5]
const _logTeff = range(convert(dtype, 3.41), convert(dtype, 4.85); step=convert(dtype, 0.01))
const _logg = range(convert(dtype, -0.5), convert(dtype, 6.0); step=convert(dtype, 0.5))

""" Unique values for dependent variables in the YBC.PHOENIX bolometric correction grid. Note that `logg = 6` is missing for [M/H] = -2.5, -3, -3.5, -4 """
const gridinfo = (logTeff = _logTeff,
                  logg = _logg,
                  MH = _mh,
                  Av = _Av,
                  Rv = dtype[3.1])
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
    return (MH = parse(dtype, mh), α_fe = parse(dtype, α_fe))
end

#########################################################

"""
    PHOENIXYBCGrid(grid::AbstractString)

Load and return the YBC PHOENIX bolometric corrections for the given photometric system `grid`,
which must be a valid entry in `BolometricCorrections.YBC.systems`.
This type is used to create instances of [`PHOENIXYBCTable`](@ref) that have fixed dependent
grid variables (\\[M/H\\], Av). This can be done either by calling an instance of
`PHOENIXYBCGrid` with `(mh, Av)` arguments or by using the appropriate constructor for [`PHOENIXYBCTable`](@ref).

```jldoctest
julia> grid = PHOENIXYBCGrid("acs_wfc")
YBC PHOENIX bolometric correction grid for photometric system YBC/acs_wfc.

julia> grid(-1.01, 0.11) # Can be called to construct table with interpolated [M/H], Av
YBC PHOENIX BT-Settl bolometric correction table with for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.11

julia> chemistry(grid) isa BolometricCorrections.MIST.MISTChemistry # Same chemical mixture as MIST
true
```
"""
struct PHOENIXYBCGrid{A <: Number, C <: AbstractVector{A}, N} <: AbstractBCGrid{A}
    data::Array{A, 5} # A should be Float32
    mag_zpt::C
    systems::Vector{String}
    name::String
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end

function PHOENIXYBCGrid(data::Array{A, 5}, mag_zpt::AbstractArray{<:Number}, systems, name, filternames) where {A}
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
    # Sort files by [M/H] value
    idxs = sortperm([_parse_filename(file).MH for file in files])
    files = files[idxs]
    # Check that MH values as parsed from filenames match gridinfo.MH
    if ~mapreduce(isapprox, &, gridinfo.MH, _parse_filename(file).MH for file in files)
        error("File [M/H] values not as expected -- please report.")
    end

    # data = zeros(dtype, length(gridinfo.logg), length(gridinfo.logTeff), length(filternames), length(files), length(gridinfo.Av))
    # for i in eachindex(files)
    #     file = files[i]
    #     FITS(file, "r") do f
    #         _logg = read(f[2], "logg")
    #         _logTeff = read(f[2], "logTeff")
    #         for j in eachindex(gridinfo.Av)
    #             Av = gridinfo.Av[j]
    #             # Figure out FITS file column name for given Av
    #             Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(Av)
    #             # data[i,j] = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
    #             tmpdata = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
    #             for idx in eachindex(_logg)
    #                 idx_lg = searchsortedfirst(gridinfo.logg, _logg[idx])
    #                 idx_lt = searchsortedfirst(gridinfo.logTeff, _logTeff[idx])
    #                 data[idx_lg, idx_lt, :, i, j] .= @view(tmpdata[idx, :])
    #             end
    #         end
    #     end
    # end

    # # Loop and fill in missing (bad) values == 0
    # for i=axes(data, 3), j=axes(data, 4), k=axes(data, 5)
    #     tmpdata = @view(data[:,:,i,j,k])
    #     if zero(dtype) in tmpdata
    #         tmpdata .= fill_bad_values(tmpdata; isbad = Base.Fix1(==, zero(dtype)), window = 1, diag = false)
    #     end
    # end

    data = zeros(dtype, length(gridinfo.logg), length(gridinfo.logTeff), length(filternames), length(files), length(gridinfo.Av))
    for i in eachindex(files)
        file = files[i]
        FITS(file, "r") do f
            _logg = read(f[2], "logg")
            _logTeff = read(f[2], "logTeff")
            u_logg = sort(unique(_logg))
            u_logTeff = sort(unique(_logTeff))
            lg_1, lg_2 = searchsortedfirst(gridinfo.logg, first(u_logg)), searchsortedfirst(gridinfo.logg, last(u_logg))
            lt_1, lt_2 = searchsortedfirst(gridinfo.logTeff, first(u_logTeff)), searchsortedfirst(gridinfo.logTeff, last(u_logTeff))
            for j in eachindex(gridinfo.Av)
                Av = gridinfo.Av[j]
                # Figure out FITS file column name for given Av
                Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(Av)
                for k in eachindex(filternames)
                    tmpdata = read(f[2], String(filternames[k])*Av_prefix)
                    data[lg_1:lg_2, lt_1:lt_2, :, i, j] .= reshape(tmpdata, length(u_logg), length(u_logTeff))
                end
            end
            # Extrapolate matrix if data does not cover full gridinfo range
            if lg_1 != 1
                for k in 1:lg_1
                    data[k, :, :, :, :] .= @view(data[lg_1, :, :, :, :])
                end
            end
            if lg_2 != lastindex(gridinfo.logg)
                for k in lg_2:lastindex(gridinfo.logg)
                    data[k, :, :, :, :] .= @view(data[lg_2, :, :, :, :])
                end
            end
            if lt_1 != 1
                for k in 1:lt_1
                    data[:, k, :, :, :] .= @view(data[:, lt_1, :, :, :])
                end
            end
            if lt_2 != lastindex(gridinfo.logTeff)
                for k in lt_2:lastindex(gridinfo.logTeff)
                    data[:, k, :, :, :] .= @view(data[:, lt_2, :, :, :])
                end
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
Base.extrema(::Type{<:PHOENIXYBCGrid}) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                          logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                          MH = (first(gridinfo.MH), last(gridinfo.MH)),
                                          Av = (first(gridinfo.Av), last(gridinfo.Av)),
                                          Rv = (first(gridinfo.Rv), last(gridinfo.Rv)))
filternames(grid::PHOENIXYBCGrid) = grid.filters
gridname(::Type{<:PHOENIXYBCGrid}) = "YBC-PHOENIX"
chemistry(::Type{<:PHOENIXYBCGrid}) = MISTChemistry()
# zeropoints(::PHOENIXYBCGrid) = zpt


#########################################################
# A single BC table, with fixed [M/H] and Av

"""
    PHOENIXYBCTable(grid::PHOENIXYBCGrid, mh::Real, Av::Real)

Interpolates the YBC PHOENIX bolometric corrections in `grid` to a fixed value of \\[M/H\\]
(`mh`) and V-band extinction (`Av`), leaving only `Teff` and `logg` as dependent
variables (the YBC PHOENIX BCs have only one `Rv` value). Returns an instance that is callable
with arguments `(Teff [K], logg [cgs])` to interpolate the bolometric corrections as a function
of temperature and surface gravity.

    PHOENIXYBCTable(grid::AbstractString, mh::Real, Av::Real)

Loads the data necessary to construct the BC table for the provided `grid` (e.g., `"acs_wfc"`)
at \\[M/H\\] = `mh` and V-band extinction `Av`. This method does not support interpolation 
in metallicity or extinction, so the arguments `mh` and `Av` must be among the values 
provided by PHOENIX (see `BolometricCorrections.YBC.PHOENIX.gridinfo.MH`).

```jldoctest
julia> grid = PHOENIXYBCGrid("acs_wfc")
YBC PHOENIX bolometric correction grid for photometric system YBC/acs_wfc.

julia> table = PHOENIXYBCTable(grid, -1.01, 0.011) # Interpolate table from full grid
YBC PHOENIX BT-Settl bolometric correction table with for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.011

julia> length(table(2755, 0.01)) == 12 # Returns BC in each filter
true

julia> size(table([2755, 2756], [0.01, 0.02])) # `table(array, array)` is also supported
(12, 2)

julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

julia> table(Table, [2755, 2756], [0.01, 0.02]) isa Table
true

julia> chemistry(table) isa BolometricCorrections.MIST.MISTChemistry # Same chemical mixture as MIST
true

julia> PHOENIXYBCTable("acs_wfc", -2.0, 0.5) isa PHOENIXYBCTable # Can construct without a PHOENIXYBCGrid
true
```
"""
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
    T = dtype # T = promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
    return PHOENIXYBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), itp, filters)
end
chemistry(::Type{<:PHOENIXYBCTable}) = MISTChemistry()
Base.show(io::IO, z::PHOENIXYBCTable) = print(io, "YBC PHOENIX BT-Settl bolometric correction table with for system $(z.name) with [M/H] ",
                                              z.MH, " and V-band extinction ", z.Av)
filternames(table::PHOENIXYBCTable) = table.filters
gridname(::Type{<:PHOENIXYBCTable}) = "YBC-PHOENIX"
# zeropoints(table::PHOENIXYBCTable) = table.mag_zpt
MH(t::PHOENIXYBCTable) = t.MH
Z(t::PHOENIXYBCTable) = Z(chemistry(t), MH(t))

# Interpolations uses `bounds` to return interpolation domain
# We will just use the hard-coded grid bounds; extremely fast
Base.extrema(::Type{<:PHOENIXYBCTable}) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                           logg = (first(gridinfo.logg), last(gridinfo.logg)))
(table::PHOENIXYBCTable)(Teff::Real, logg::Real) = table(promote(Teff, logg)...)
(table::PHOENIXYBCTable)(Teff::T, logg::T) where {T <: Real} = table.itp(logg, log10(Teff))
# Data are naturally Float32 -- convert hardware numeric args for faster evaluation and guarantee Float32 output
(table::PHOENIXYBCTable)(Teff::HardwareNumeric, logg::HardwareNumeric) = table(convert(dtype, Teff), convert(dtype, logg))
# to broadcast over both teff and logg, you do table.(teff, logg')

function PHOENIXYBCTable(grid::AbstractString, mh::Real, Av::Real; prefix::AbstractString="YBC")
    grid, prefix = String(grid), String(prefix)
    check_prefix(prefix)
    @argcheck mapreduce(isapprox(mh), |, gridinfo.MH) "Provided [M/H] $mh not in available values $(gridinfo.MH); use PHOENIXYBCGrid for grid interpolation."
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use PHOENIXYBCGrid for grid interpolation."
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
        # println(length(unique(read(f[2], "logg"))) * length(unique(read(f[2], "logTeff"))))
        logg = mh > -2.5 ? gridinfo.logg : gridinfo.logg[begin:end-1]
        newdata = repack_submatrix(data, length(logg), length(gridinfo.logTeff), Val(length(filternames)))
        itp = cubic_spline_interpolation((logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
        return PHOENIXYBCTable(mh, Av, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
    end
end

function PHOENIXYBCTable(grid::PHOENIXYBCGrid, mh::Real, Av::Real)
    check_vals(mh, Av, gridinfo)
    filters = filternames(grid)
    data = grid.data

    Av_vec = SVector{length(gridinfo.Av), dtype}(gridinfo.Av) # Need vector to use searchsortedfirst
    MH_vec = gridinfo.MH # SVector(gridinfo.MH)

    # Exact values are in grid; no interpolation necessary
    if mh ∈ gridinfo.MH && Av ∈ gridinfo.Av
        submatrix = data[:, :, :, searchsortedfirst(MH_vec, mh), searchsortedfirst(Av_vec, Av)]
    else
        if mh ∈ gridinfo.MH
            MH_idx = searchsortedfirst(MH_vec, mh)
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            mat1 = data[:, :, :, MH_idx, Av_idx]
            mat2 = data[:, :, :, MH_idx, Av_idx + 1]
            submatrix = interp1d(Av, Av_vec[Av_idx], Av_vec[Av_idx + 1], mat1, mat2)
        elseif Av ∈ gridinfo.Av
            Av_idx = searchsortedfirst(Av_vec, Av)
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mat1 = data[:, :, :, MH_idx, Av_idx]
            mat2 = data[:, :, :, MH_idx + 1, Av_idx]
            submatrix = interp1d(mh, MH_vec[MH_idx], MH_vec[MH_idx + 1], mat1, mat2)
        else
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            Av1, Av2 = Av_vec[Av_idx], Av_vec[Av_idx + 1]
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mh1, mh2 = MH_vec[MH_idx], MH_vec[MH_idx + 1]

            mat1_1 = data[:, :, :, MH_idx, Av_idx]
            mat2_1 = data[:, :, :, MH_idx + 1, Av_idx]
            mat1_2 = data[:, :, :, MH_idx, Av_idx + 1]
            mat2_2 = data[:, :, :, MH_idx + 1, Av_idx + 1]
            # Perform bilinear interpolation
            submatrix = interp2d(mh, Av, mh1, mh2, Av1, Av2, mat1_1, mat2_1, mat1_2, mat2_2)
        end
    end
    newdata = repack_submatrix(submatrix, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filters)))
    itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
    return PHOENIXYBCTable(mh, Av, grid.mag_zpt, grid.systems, grid.name, itp, filters)
end
PHOENIXYBCTable(grid::PHOENIXYBCGrid, mh::HardwareNumeric, Av::HardwareNumeric) = PHOENIXYBCTable(grid, convert(dtype, mh), convert(dtype, Av))
end # module