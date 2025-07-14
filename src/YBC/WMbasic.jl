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
using ..YBC: dtype, pull_table, parse_filterinfo, check_prefix, check_vals

# export ATLAS9YBCTable, ATLAS9YBCGrid, ATLAS9Chemistry

""" `NTuple{6, Symbol}` listing the dependent variables in the YBC.WMbasic BC grid. """
const _dependents = (:logTeff, :logg, :Z, :Mdot, :Av, :Rv)
""" A_v values in files. For each filter "J", each fits file will have columns "J", "J_Av0.5", "J_Av1", and so on."""
const _Av = (0, 0.5, 1, 2, 5, 10, 20) # Mix of float and integer makes parsing FITS columns easier later
""" Unique values of metal mass fraction Z for the ATLAS9 models. """
const _Z = dtype[0.0001, 0.0005, 0.001, 0.004, 0.008, 0.02]
"""Unique values of mass loss range in solar masses per year."""
const _Mdot = logrange(convert(dtype, 1e-7), convert(dtype, 1e-5), 3) # dtype[1e-5, 1e-6, 1e-7]
const _logTeff = range(convert(dtype, 4.325), convert(dtype, 5.0); step=convert(dtype, 0.025))
const _logg = range(convert(dtype, 2.5), convert(dtype, 5.0); step=convert(dtype, 0.5))

"""Unique values for dependent variables in the YBC.ATLAS9 bolometric correction grid."""
const gridinfo = (logTeff = _logTeff,
                  logg = _logg,
                  Z = _Z,
                  # MH = _mh,
                  Mdot = _Mdot,
                  Av = _Av,
                  Rv = dtype[3.1])
@compat public gridinfo

"""
    _parse_filename(f::AbstractString)

Return Z and Mdot of YBC WMBasic model given a filename (example: "Avodonnell94Rv3.1WM_Z0.0001Mdot5.BC.fits").
"""
function _parse_filename(f::AbstractString)
    f = basename(f)
    @argcheck occursin("WM_", f)
    Z = parse(dtype, split(split(f, "WM_Z")[2], "Mdot")[1])
    Mdot = exp10(-1 * parse(dtype, split(split(f, "Mdot")[2], ".BC")[1]))
    return (Z = Z, Mdot = Mdot)
end


function WMbasicYBCTable(grid::AbstractString, zval::Real, Av::Real, Mdot::Real; prefix::AbstractString="YBC")
    grid, prefix = String(grid), String(prefix)
    check_prefix(prefix)
    @argcheck mapreduce(isapprox(zval), |, gridinfo.Z) "Provided metal mass fraction $zval not in available values $(gridinfo.Z); use WMbasicYBCGrid for grid interpolation."
    @argcheck mapreduce(isapprox(Av), |, gridinfo.Av) "Provided Av $Av not in available values $(gridinfo.Av); use WMbasicYBCGrid for grid interpolation."
    path = pull_table(grid, prefix)
    files = filter(x->occursin("WM_", x), readdir(joinpath(path, "regrid"); join=true))
    if length(files) == 0
        error("""No files found for grid $grid in the given YBC directory $prefix. prefix="YBC" has the greatest number of filters and is recommended.""")
    elseif length(files) != length(gridinfo.Z) * length(gridinfo.Mdot)
        error("Number of files found for grid $grid is $(length(files)), expected $(length(gridinfo.Z)). Data may be corrupted. \\
        Recommend purging data with `BolometricCorrections.YBC.remove_table($grid; prefix = $prefix)` and rerunning.")
    end
    filterinfo = parse_filterinfo(joinpath(path, "filter.info"))
    filternames = filterinfo.names

    # Figure out which file we need for given [M/H]
    goodfile = files[findfirst(x -> (x.Z ≈ zval) & (x.Mdot ≈ Mdot), _parse_filename(file) for file in files)]
    # Figure out FITS file column name for given Av
    Av_prefix = isapprox(0, Av) ? "" : "_Av" * string(gridinfo.Av[findfirst(≈(Av), gridinfo.Av)])
    # Need to loop over FITS files for different Mdot, write into extra dimension of data
    data = Array{dtype, 3}(undef, length(gridinfo.logg), length(gridinfo.logTeff), length(gridinfo.Mdot))
    for i in eachindex(gridinfo.Mdot)
        mdot_str = @sprintf("%1i", log10(gridinfo.Mdot[i]))[2] # Returns something like "-5" for 1f-5, only want "5"
        return mdot_str
    # FITS(goodfile, "r") do f
    #     data = reduce(hcat, read(f[2], String(filt)*Av_prefix) for filt in filternames)
    #     # Pack data into (length(logg), length(logTeff)) Matrix{SVector} for interpolation
    #     newdata = repack_submatrix(data, length(gridinfo.logg), length(gridinfo.logTeff), Val(length(filternames)))
    #     itp = cubic_spline_interpolation((gridinfo.logg, gridinfo.logTeff), newdata; extrapolation_bc=Flat())
    #     return itp
    #     return WMbasicYBCTable(zval, Av, Mdot, filterinfo.mag_zeropoint, String.(filterinfo.photometric_system), prefix*"/"*grid, itp, tuple(Symbol.(filternames)...))
    # end
    end
end


end # module