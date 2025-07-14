"""
YBC submodule exposing the bolometric corrections computed by the YBC team based on their own WM-basic models. These models assume the PARSEC solar abundances [Bressan2012](@citet).
"""
module WMbasic

using ArgCheck: @argcheck
using Compat: @compat, logrange
using FITSIO: FITS
using Interpolations: cubic_spline_interpolation, Throw, Flat
using Printf: @sprintf
using StaticArrays: SVector


using ...BolometricCorrections: repack_submatrix, AbstractBCTable, AbstractBCGrid, interp1d, interp2d
import ...BolometricCorrections: zeropoints, filternames, chemistry # , Y_p, X, X_phot, Y, Y_phot, Z, Z_phot, MH # vegamags, abmags, stmags, Mbol, Lbol
using ..YBC: dtype, pull_table, parse_filterinfo, check_prefix, check_vals, PARSECChemistry, Z, MH

# export ATLAS9YBCTable, ATLAS9YBCGrid, ATLAS9Chemistry

""" `NTuple{6, Symbol}` listing the dependent variables in the YBC.WMbasic BC grid. """
const _dependents = (:logTeff, :logg, :Z, :Mdot, :Av, :Rv)
""" A_v values in files. For each filter "J", each fits file will have columns "J", "J_Av0.5", "J_Av1", and so on."""
const _Av = (0, 0.5, 1, 2, 5, 10, 20) # Mix of float and integer makes parsing FITS columns easier later
""" Unique values of metal mass fraction Z for the ATLAS9 models. """
const _Z = dtype[0.0001, 0.0005, 0.001, 0.004, 0.008, 0.02]
"""Unique values of mass loss range in solar masses per year."""
const _Mdot = logrange(convert(dtype, 1e-7), convert(dtype, 1e-5), 3) # dtype[1e-5, 1e-6, 1e-7]
const _logTeff = range(convert(dtype, 4.3), convert(dtype, 5.0); step=convert(dtype, 0.025))
const _logg = range(convert(dtype, 2.5), convert(dtype, 6.0); step=convert(dtype, 0.5))
# Mdot = 1.0e-5 only goes to logg = 5.5, not 6.0
# Z = 0.0001 has logTeff = 4.3:0.025:5, logg = 2.5:0.5:5.5 except for Mdot=1e-5 which has logg = 2.5:0.5:5, logTeff = 4.325:0.025:5
# Z = 0.0005 has logTeff = 4.3:0.025:5, logg = 2.5:0.5:6.0 for Mdot=1e-7, logg = 2.5:0.5:5.5 for Mdot=1e-6, logg = 2.5:0.5:5.0 for Mdot=1e-5
# Z = 0.001 same as above
# Z = 0.004 same as above, but logg = 2.5:0.5:5.5 for Mdot=1e-5
# Z = 0.008 logTeff as above, logg=3:0.5:6 for Mdot=1e-7, logg=2.5:0.5:5.5 for Mdot=1e-6, 1e-5
# Z = 0.02 logTeff as above, logg=2.5:0.5:6.0 for Mdot=1e-7 and 1e-6, logg=2.5:0.5:5.5 for Mdot=1e-5

"""Unique values for dependent variables in the YBC.ATLAS9 bolometric correction grid."""
const gridinfo = (logTeff = _logTeff,
                  logg = _logg,
                  Z = _Z,
                  MH = [convert(dtype, MH(PARSECChemistry(), z)) for z in _Z],
                  Mdot = _Mdot,
                  Av = _Av,
                  Rv = dtype[3.1])
@compat public gridinfo

"""
    _parse_filename(f::AbstractString)

Return Z and Mdot of YBC WMbasic model given a filename (example: "Avodonnell94Rv3.1WM_Z0.0001Mdot5.BC.fits").
"""
function _parse_filename(f::AbstractString)
    f = basename(f)
    @argcheck occursin("WM_", f)
    Z = parse(dtype, split(split(f, "WM_Z")[2], "Mdot")[1])
    Mdot = exp10(-1 * parse(dtype, split(split(f, "Mdot")[2], ".BC")[1]))
    return (Z = Z, Mdot = Mdot)
end

#########################################################

"""
    WMbasicYBCGrid(grid::AbstractString)

Load and return the YBC WM-basic bolometric corrections for the given photometric system `grid`,
which must be a valid entry in `BolometricCorrections.YBC.systems`.
This type is used to create instances of [`WMbasicYBCTable`](@ref) that have fixed dependent
grid variables (\\[M/H\\], Av). This can be done either by calling an instance of
`WMbasicYBCGrid` with `(mh, Av)` arguments or by using the appropriate constructor for [`WMbasicYBCTable`](@ref).

```jldoctest
julia> grid = WMbasicYBCTable("acs_wfc")
YBC WM-basic bolometric correction grid for photometric system YBC/acs_wfc.

julia> grid(-1.01, 0.11) # Can be called to construct table with interpolated [M/H], Av
YBC WM-basic bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.11
```
"""
struct WMbasicYBCGrid{A <: Number, C <: AbstractVector{A}, N} <: AbstractBCGrid{A}
    data::Array{Matrix{A}, 3} # A should be Float32
    mag_zpt::C
    systems::Vector{String}
    name::String
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end

function WMbasicYBCGrid(data::Array{Matrix{A}, 3}, mag_zpt::AbstractArray{<:Number}, systems, name, filternames) where {A}
    return WMbasicYBCGrid(data, convert.(A, mag_zpt), String.(systems), String(name), tuple(Symbol.(filternames)...))
end

function WMbasicYBCGrid(grid::AbstractString; prefix::AbstractString="YBC")
    check_prefix(prefix)
    path = pull_table(String(grid), String(prefix))
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Read all data and pack into dense matrix
    data = Array{Matrix{dtype}, 3}(undef, length(gridinfo.Z), length(gridinfo.Mdot), length(gridinfo.Av))
    for i in eachindex(gridinfo.Z)
        z = gridinfo.Z[i]
        # println("-----------------------")
        # println("Z: ", z)
        for j in eachindex(gridinfo.Mdot)
            mdot = gridinfo.Mdot[j]
            println("Z: ", z, " Mdot: ", mdot)
            file = joinpath(path, "regrid", "Avodonnell94Rv3.1WM_Z" * string(z) * "Mdot" * @sprintf("%1i", log10(mdot))[2] * ".BC.fits")
            if !isfile(file)
                error("YBC WM-basic file $file missing. Data may be corrupted. Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
            end
            FITS(file, "r") do f
                for k in eachindex(gridinfo.Av)
                    Av = gridinfo.Av[k]
                    # println("Av: ", Av)
                    # Figure out FITS file column name for given Av
                    Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(Av)
                    logg = read(f[2], "logg")
                    logTeff = read(f[2], "logTeff")

                    # mat of correct size that we will save
                    # mat1 = Matrix{dtype}(undef, length(gridinfo.logg), length(gridinfo.logTeff))
                    mat1 = zeros(dtype, length(gridinfo.logg) * length(gridinfo.logTeff), length(filternames))
                    # mat of irregular size due to irregular sampling of logteff, logg
                    mat2 = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)

                    for ii in eachindex(logg, logTeff)
                        mati = searchsortedfirst(gridinfo.logg, logg[ii])
                        matj = searchsortedfirst(gridinfo.logTeff, logTeff[ii])
                        mat1[mati + ((matj - 1) * length(gridinfo.logg)), :] .= mat2[ii, :]
                    end

                    # 4 special cases; logg=2.5 missing, logg=5.5 missing, logg=6.0 missing, logTeff = 4.3 missing
                    u_logg = unique(logg)
                    u_logTeff = unique(logTeff)
                    # for lgi in findall(lg ∉ u_logg for lg in gridinfo.logg)
                    #     # println([lgi + length(gridinfo.logg) * N for N in 0:(length(gridinfo.logTeff))])
                    #     for ii in [lgi + length(gridinfo.logg) * N for N in 0:(length(gridinfo.logTeff) - 1)]
                    #         # println(ii)
                    #         if checkbounds(Bool, mat1, ii + 1, 1)
                    #             # Use same temperature, next value of logg (for logg=2.5 missing)
                    #             # println("yep1")
                    #             jj = findfirst(!=(0), mat1[ii:end, 1])
                    #             if isnothing(jj)
                    #                 jj = ii - 1
                    #             else
                    #                 jj += ii
                    #             end
                    #             mat1[ii, :] .= mat1[jj, :]
                    #             # if mat1[ii + 1, 1] != 0
                    #             #     mat1[ii, :] .= mat1[ii + 1, :]
                    #             # else
                    #             #     mat1[ii, :] .= mat1[ii - 1, :]
                    #             # end
                    #         elseif checkbounds(Bool, mat1, ii - 1, 1)
                    #             # Use same temperature, last value of logg (for log=5.5, 6.0 missing)
                    #             println("yep2")
                    #             mat1[ii, :] .= mat1[ii - 1, :]
                    #         else
                    #             println("failed")
                    #         end
                    #     end
                    # end
                    if 2.5f0 ∉ u_logg
                        for ii in [1 + length(gridinfo.logg) * N for N in 0:(length(gridinfo.logTeff) - 1)]
                            # Average of same temperature, next logg and same logg, next temperature
                            # mat1[ii, :] .= (mat1[2, :] .+ mat1[ii + length(gridinfo.logg), :]) ./ 2
                            # Actually, just use same temperature, next logg
                            mat1[ii, :] .= mat1[ii + 1, :]
                        end
                    end
                    if 5.5f0 ∉ u_logg
                        for ii in [7 + length(gridinfo.logg) * N for N in 0:(length(gridinfo.logTeff) - 1)]
                            mat1[ii, :] .= mat1[ii - 1, :] # Use same temperature, last logg
                        end
                    end
                    if 6.0f0 ∉ u_logg
                        for ii in [8 + length(gridinfo.logg) * N for N in 0:(length(gridinfo.logTeff) - 1)]
                            mat1[ii, :] .= mat1[ii - 1, :] # Use same temperature, last logg
                        end
                    end
                    if 4.3f0 ∉ u_logTeff
                        # for ii in [1 + length(gridinfo.logTeff) * N for N in 0:(length(gridinfo.logg) - 1)]
                        for ii in 1:length(gridinfo.logg)
                            mat1[ii, :] .= mat1[ii + length(gridinfo.logg), :] # Use next temperature, same logg
                        end
                    end

                    # println(size(mat))
                    # println("logg: ", unique(read(f[2], "logg")))
                    # println("logTeff: ", unique(read(f[2], "logTeff")))
                    data[i,j,k] = mat1
                end
            end
        end
    end
    return WMbasicYBCGrid(data, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, filternames)
end
(grid::WMbasicYBCGrid)(mh::Real, Av::Real) = WMbasicYBCTable(grid, mh, Av)
Base.show(io::IO, z::WMbasicYBCGrid) = print(io, "YBC WM-basic bolometric correction grid for photometric system $(z.name).")
# # function Table(grid::WMbasicYBCGrid)
# #     data = grid.data
# #     tables = Vector{Table}(undef, length(data))
# # end
Base.extrema(::WMbasicYBCGrid) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                  logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                  # MH = (first(gridinfo.MH), last(gridinfo.MH)),
                                  Z = (first(gridinfo.Z), last(gridinfo.Z)), 
                                  Av = (first(gridinfo.Av), last(gridinfo.Av)),
                                  Mdot = (first(gridinfo.Mdot), last(gridinfo.Mdot)),
                                  Rv = (first(gridinfo.Rv), last(gridinfo.Rv)))
filternames(grid::WMbasicYBCGrid) = grid.filters
# zeropoints(::WMbasicYBCGrid) = zpt
chemistry(::Type{WMbasicYBCGrid}) = PARSECChemistry()

#########################################################
# A single BC table, with fixed [M/H] and Av

"""
    WMbasicYBCTable(grid::WMbasicYBCGrid, mh::Real, Av::Real)

Interpolates the YBC WM-basic bolometric corrections in `grid` to a fixed value of \\[M/H\\]
(`mh`), V-band extinction (`Av`), leaving only `Teff`, `logg`, and `Mdot as dependent
variables (the YBC WM-basic BCs have only one `Rv` value).
Returns an instance that is callable with arguments `(Teff [K], logg [cgs], 
Mdot [solMass / yr])` to interpolate the bolometric corrections as a function
of temperature, surface gravity, and outflow rate.

```jldoctest
julia> grid = WMbasicYBCGrid("acs_wfc")
YBC WM-basic bolometric correction grid for photometric system YBC/acs_wfc.

julia> table = WMbasicYBCGrid(grid, -1.01, 0.011) # Interpolate table from full grid
YBC WM-basic bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.011

julia> length(table(25_0254.0, 2.54, 5e-6)) == 12 # Returns BC in each filter
true

julia> size(table([25_0254.0, 25_0354.0], [2.54, 2.56], [4e-6, 5e-6])) # `table(array, array)` is also supported
(12, 2)

julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

julia> table(Table, [25_0254.0, 25_0354.0], [2.54, 2.56], [4e-6, 5e-6]) isa Table
true
```
"""
struct WMbasicYBCTable{A <: Real, B, N} <: AbstractBCTable{A}
    MH::A
    Av::A
    mag_zpt::Vector{A}
    systems::Vector{String}
    name::String
    itp::B     # Interpolator object
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
function WMbasicYBCTable(MH::Real, Av::Real, mag_zpt::Vector{<:Real}, systems, name, itp, filters)
    T = dtype # promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
    return WMbasicYBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), itp, filters)
end
Base.show(io::IO, z::WMbasicYBCTable) = print(io, "YBC WM-basic bolometric correction table with for system $(z.name) with [M/H] ",
                                              z.MH, " and V-band extinction ", z.Av)
filternames(table::WMbasicYBCTable) = table.filters
# zeropoints(table::WMbasicYBCTable) = table.mag_zpt
chemistry(::Type{WMbasicYBCTable}) = PARSECChemistry()

# Interpolations uses `bounds` to return interpolation domain
# We will just use the hard-coded grid bounds; extremely fast
Base.extrema(::WMbasicYBCTable) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
                                   logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                   Mdot = (first(gridinfo.Mdot), last(gridinfo.Mdot)))
(table::WMbasicYBCTable)(Teff::Real, logg::Real, Mdot::Real) = table.itp(logg, log10(Teff), Mdot)
# Data are naturally Float32 -- convert Float64 args for faster evaluation (~35% faster)
(table::WMbasicYBCTable)(Teff::Float64, logg::Float64, Mdot::Float64) = table(convert(dtype, Teff), convert(dtype, logg), convert(dtype, Mdot))
# to broadcast over both teff and logg, you do table.(teff, logg')

function WMbasicYBCTable(grid::WMbasicYBCGrid, mh::Real, Av::Real)
    check_vals(mh, Av, gridinfo)
    filters = filternames(grid)
    data = grid.data

    Av_vec = SVector(gridinfo.Av) # Need vector to use searchsortedfirst
    MH_vec = gridinfo.MH

    if mh ∈ gridinfo.MH && Av ∈ gridinfo.Av
        # Exact values are in grid; no interpolation necessary
        submatrix = data[searchsortedfirst(MH_vec, mh), :, searchsortedfirst(Av_vec, Av)]
    else
        if mh ∈ gridinfo.MH
            MH_idx = searchsortedfirst(MH_vec, mh)
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            vecmat1 = data[MH_idx, :, Av_idx] # Vector{Matrix} of length gridinfo.Mdot
            vecmat2 = data[MH_idx, :, Av_idx + 1]
            submatrix = [interp1d(Av, Av_vec[Av_idx], Av_vec[Av_idx + 1], vecmat1[i], vecmat2[i]) for i in eachindex(vecmat1, vecmat2)]
        elseif Av ∈ gridinfo.Av
            Av_idx = searchsortedfirst(Av_vec, Av)
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            vecmat1 = data[MH_idx, :, Av_idx]
            vecmat2 = data[MH_idx + 1, :, Av_idx]
            submatrix = [interp1d(mh, MH_vec[MH_idx], MH_vec[MH_idx + 1], vecmat1[i], vecmat2[i]) for i in eachindex(vecmat1, vecmat2)]
        else
            Av_idx = searchsortedfirst(Av_vec, Av) - 1
            Av1, Av2 = Av_vec[Av_idx], Av_vec[Av_idx + 1]
            MH_idx = searchsortedfirst(MH_vec, mh) - 1
            mh1, mh2 = MH_vec[MH_idx], MH_vec[MH_idx + 1]

            vecmat1_1 = data[MH_idx, :, Av_idx]
            vecmat2_1 = data[MH_idx + 1, :, Av_idx]
            vecmat1_2 = data[MH_idx, :, Av_idx + 1]
            vecmat2_2 = data[MH_idx + 1, :, Av_idx + 1]
            # Perform bilinear interpolation
            submatrix = [interp2d(mh, Av, mh1, mh2, Av1, Av2, vecmat1_1[i], vecmat2_1[i], vecmat1_2[i], vecmat2_2[i]) for i in eachindex(vecmat1_1)]
        end
    end
    return submatrix
    # Need to check dimensionality for repack_submatrix
    # newdata = repack_submatrix(submatrix, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filters)))
    # itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
    # return ATLAS9YBCTable(mh, Av, grid.mag_zpt, grid.systems, grid.name, itp, filters)
end
WMbasicYBCTable(grid::WMbasicYBCGrid, mh::Float64, Av::Float64) = WMbasicYBCTable(grid, convert(dtype, mh), convert(dtype, Av))




# function WMbasicYBCTable(grid::AbstractString, zval::Real, Av::Real, Mdot::Real; prefix::AbstractString="YBC")
#     grid, prefix = String(grid), String(prefix)
#     check_prefix(prefix)
#     @argcheck mapreduce(isapprox(zval), |, gridinfo.Z) "Provided metal mass fraction $zval not in available values $(gridinfo.Z); use WMbasicYBCGrid for grid interpolation."
#     @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use WMbasicYBCGrid for grid interpolation."
#     path = pull_table(grid, prefix)
#     files = filter(x->occursin("WM_", x), readdir(joinpath(path, "regrid"); join=true))
#     if length(files) == 0
#         error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
#     elseif length(files) != length(gridinfo.Z) * length(gridinfo.Mdot)
#         error("Number of files found for grid $grid is $(length(files)), expected $(length(gridinfo.Z)). Data may be corrupted. \\
#         Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
#     end
#     filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
#     filternames = filterinfo.names

#     # Figure out which file we need for given [M/H]
#     goodfile = files[findfirst(x -> (x.Z ≈ zval) & (x.Mdot ≈ Mdot), _parse_filename(file) for file in files)]
#     # Figure out FITS file column name for given Av
#     Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
#     FITS(goodfile, "r") do f
#         data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
#         # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
#         newdata = repack_submatrix(data, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filternames)))
#         itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
#         return WMbasicYBCTable(zval, Av, Mdot, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
#     end
# end

    # # Need to loop over FITS files for different Mdot, write into extra dimension of data
    # data = Array{dtype, 3}(undef, length(gridinfo.logg), length(gridinfo.logTeff), length(gridinfo.Mdot))
    # for i in eachindex(gridinfo.Mdot)
    #     mdot_str = @sprintf("%1i", log10(gridinfo.Mdot[i]))[2] # Returns something like "-5" for 1f-5, only want "5"
    # end

#########################################################
# A single BC table, with fixed [M/H] and Av

# """
#     WMbasicYBCTable(grid::ATLAS9YBCGrid, mh::Real, Av::Real)

# Interpolates the YBC ATLAS9 bolometric corrections in `grid` to a fixed value of \\[M/H\\]
# (`mh`), V-band extinction (`Av`), leaving only `Teff`, `logg`, and `Mdot as dependent
# variables (the YBC WM-basic BCs have only one `Rv` value). 
# Returns an instance that is callable with arguments `(Teff [K], logg [cgs], 
# Mdot [solMass / yr])` to interpolate the bolometric corrections as a function
# of temperature, surface gravity, and outflow rate.

# ```jldoctest
# julia> grid = WMbasicYBCGrid("acs_wfc")
# YBC WM-basic bolometric correction grid for photometric system YBC/acs_wfc.

# julia> table = WMbasicYBCGrid(grid, -1.01, 0.011) # Interpolate table from full grid
# YBC WM-basic bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.011

# julia> length(table(25_0254.0, 2.54)) == 12 # Returns BC in each filter
# true

# julia> size(table([25_0254.0, 25_0354.0], [2.54, 2.56])) # `table(array, array)` is also supported
# (12, 2)

# julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

# julia> table(Table, [25_0254.0, 25_0354.0], [2.54, 2.56]) isa Table
# true
# ```
# """
# struct ATLAS9YBCTable{A <: Real, B, N} <: AbstractBCTable{A}
#     MH::A
#     Av::A
#     mag_zpt::Vector{A}
#     systems::Vector{String}
#     name::String
#     itp::B     # Interpolator object
#     filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
# end
# function ATLAS9YBCTable(MH::Real, Av::Real, mag_zpt::Vector{<:Real}, systems, name, itp, filters)
#     T = dtype # promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
#     return ATLAS9YBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), itp, filters)
# end
# Base.show(io::IO, z::ATLAS9YBCTable) = print(io, "YBC ATLAS9 bolometric correction table with for system $(z.name) with [M/H] ",
#                                               z.MH, " and V-band extinction ", z.Av)
# filternames(table::ATLAS9YBCTable) = table.filters
# # zeropoints(table::ATLAS9YBCTable) = table.mag_zpt

# # Interpolations uses `bounds` to return interpolation domain
# # We will just use the hard-coded grid bounds; extremely fast
# Base.extrema(::ATLAS9YBCTable) = (Teff = (exp10(first(gridinfo.logTeff)), exp10(last(gridinfo.logTeff))), 
#                                   logg = (first(gridinfo.logg), last(gridinfo.logg)))
# (table::ATLAS9YBCTable)(Teff::Real, logg::Real) = table.itp(logg, log10(Teff))
# # Data are naturally Float32 -- convert Float64 args for faster evaluation (~35% faster)
# (table::ATLAS9YBCTable)(Teff::Float64, logg::Float64) = table(convert(dtype, Teff), convert(dtype, logg))
# # to broadcast over both teff and logg, you do table.(teff, logg')

# function WMbasicYBCTable(grid::AbstractString, zval::Real, Av::Real, Mdot::Real; prefix::AbstractString="YBC")
#     grid, prefix = String(grid), String(prefix)
#     check_prefix(prefix)
#     @argcheck mapreduce(isapprox(zval), |, gridinfo.Z) "Provided metal mass fraction $zval not in available values $(gridinfo.Z); use WMbasicYBCGrid for grid interpolation."
#     @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use WMbasicYBCGrid for grid interpolation."
#     path = pull_table(grid, prefix)
#     files = filter(x->occursin("WM_", x), readdir(joinpath(path, "regrid"); join=true))
#     if length(files) == 0
#         error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
#     elseif length(files) != length(gridinfo.Z) * length(gridinfo.Mdot)
#         error("Number of files found for grid $grid is $(length(files)), expected $(length(gridinfo.Z)). Data may be corrupted. \\
#         Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
#     end
#     filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
#     filternames = filterinfo.names

#     # Figure out which file we need for given [M/H]
#     goodfile = files[findfirst(x -> (x.Z ≈ zval) & (x.Mdot ≈ Mdot), _parse_filename(file) for file in files)]
#     # Figure out FITS file column name for given Av
#     Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
#     FITS(goodfile, "r") do f
#         data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
#         # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
#         newdata = repack_submatrix(data, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filternames)))
#         itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
#         return WMbasicYBCTable(zval, Av, Mdot, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
#     end
# end

#     # # Need to loop over FITS files for different Mdot, write into extra dimension of data
#     # data = Array{dtype, 3}(undef, length(gridinfo.logg), length(gridinfo.logTeff), length(gridinfo.Mdot))
#     # for i in eachindex(gridinfo.Mdot)
#     #     mdot_str = @sprintf("%1i", log10(gridinfo.Mdot[i]))[2] # Returns something like "-5" for 1f-5, only want "5"
#     # end

end # module