# Here we adapt the framework of MIST.jl to the new MIST v2.5 BC tables, which add [α/Fe] as a dependent variable.

######################################
# MIST v2.5 grid constants
######################################

"""
Unique values of `lgTef` (log₁₀ Teff) in the MIST v2.5 BC tables.
`lgTef` is the first column in the raw v2.5 files; use `10 .^ _mist_v2_lgTef` for linear Teff.
"""
const _mist_v2_lgTef = (3.17609, 3.17634, 3.1775, 3.17998, 3.18407, 3.19003, 3.19808, 3.20841,
                        3.22122, 3.23667, 3.25493, 3.27614, 3.30046, 3.32801, 3.35893, 3.39335,
                        3.43139, 3.47317, 3.5188, 3.5684, 3.62208, 3.67993, 3.74207, 3.8086,
                        3.87961, 3.9552, 4.03546, 4.12049, 4.21038, 4.30521, 4.40508, 4.51008,
                        4.62027, 4.73576, 4.85661, 4.98292, 5.11476, 5.2522, 5.39534, 5.54423,
                        5.69897) # length=41
""" Unique values of `logg` in the MIST v2.5 BC tables. """
const _mist_v2_logg = (-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                       5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5) # length=22
""" Unique values of [Fe/H] in the MIST v2.5 BC tables. """
const _mist_v2_feh = (-3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5,
                      -0.25, 0.0, 0.25, 0.5) # length=15
""" Unique values of [α/Fe] in the MIST v2.5 BC tables. """
const _mist_v2_afe = (-0.2, 0.0, 0.2, 0.4, 0.6) # length=5
""" Unique values of ``A_V`` in the MIST v2.5 BC tables. """
const _mist_v2_Av = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0) # length=13
""" Unique values of ``R_V`` in the MIST v2.5 BC tables. """
const _mist_v2_Rv = (3.1,) # length=1

""" Unique values for dependent variables in the MIST v2.5 bolometric correction grid. """
const gridinfov2 = (lgTef = _mist_v2_lgTef,
                    logg  = _mist_v2_logg,
                    feh   = _mist_v2_feh,
                    afe   = _mist_v2_afe,
                    Av    = _mist_v2_Av,
                    Rv    = _mist_v2_Rv)
@compat public gridinfov2

######################################
# MIST v2.5 validation and table I/O
######################################

"""
    check_vals_v2(feh, afe, Av)

Validate arguments for the MIST v2.5 BC grid. Throws `ArgumentError` on out-of-range values.
"""
function check_vals_v2(feh, afe, Av)
    feh_ext = extrema(_mist_v2_feh)
    if feh < first(feh_ext) || feh > last(feh_ext)
        throw(ArgumentError("Provided [Fe/H] $feh is outside the bounds for the MIST v2.5 BC tables $feh_ext"))
    end
    afe_ext = extrema(_mist_v2_afe)
    if afe < first(afe_ext) || afe > last(afe_ext)
        throw(ArgumentError("Provided [α/Fe] $afe is outside the bounds for the MIST v2.5 BC tables $afe_ext"))
    end
    Av_ext = extrema(_mist_v2_Av)
    if Av < first(Av_ext) || Av > last(Av_ext)
        throw(ArgumentError("Provided A_v $Av is outside the bounds for the MIST v2.5 BC tables $Av_ext"))
    end
end

# Data layout within the processed CSV:
#   outer → inner: feh, afe, Av, lgTef, logg (fastest)
# Rows per (lgTef) slice: length(_mist_v2_logg) = 22
# Rows per Av block:      length(_mist_v2_lgTef) * length(_mist_v2_logg) = 41 * 22 = 902
# Rows per (feh, afe):    length(_mist_v2_Av) * 902 = 13 * 902 = 11726
# Rows per feh:           length(_mist_v2_afe) * 11726 = 5 * 11726 = 58630
function select_subtable_v2(table::Table, feh::Real, afe::Real, Av::Real)
    @argcheck feh ∈ _mist_v2_feh && afe ∈ _mist_v2_afe && Av ∈ _mist_v2_Av
    nrows_per_Av      = length(_mist_v2_lgTef) * length(_mist_v2_logg)  # 902
    nrows_per_feh_afe = length(_mist_v2_Av)    * nrows_per_Av           # 11726
    nrows_per_feh     = length(_mist_v2_afe)   * nrows_per_feh_afe      # 58630
    idx_feh = 1 + nrows_per_feh * (findfirst(==(feh), _mist_v2_feh) - 1)
    idx_afe = idx_feh + nrows_per_feh_afe * (findfirst(==(afe), _mist_v2_afe) - 1)
    idx_Av  = idx_afe + nrows_per_Av      * (findfirst(==(Av),  _mist_v2_Av)  - 1)
    return table[idx_Av:(idx_Av + nrows_per_Av - 1)]
end

######################################
# MISTBCGridv2
######################################

"""
    MISTBCGridv2(grid::AbstractString)

Load and return the MIST v2.5 bolometric corrections for the given photometric system `grid`.
This type is used to create instances of [`MISTBCTablev2`](@ref) with fixed dependent grid
variables (\\[Fe/H\\], \\[α/Fe\\], Av). Call an instance with `(feh, afe, Av)` arguments or use
the [`MISTBCTablev2`](@ref) constructor directly.

```jldoctest
julia> grid = MISTBCGridv2("JWST")
MIST v2.5 bolometric correction grid for photometric system MIST_JWST_v2.5

julia> grid(-1.01, 0.2, 0.11)
MIST v2.5 bolometric correction table with [Fe/H] -1.01, [α/Fe] 0.2, and V-band extinction 0.11
```
"""
struct MISTBCGridv2{A, B, C <: AbstractString} <: AbstractBCGrid{A}
    table::B
    filename::C
    function MISTBCGridv2(table::B, filename::C) where {B, C}
        A = Base.promote_eltype(first(table))
        new{A, B, C}(table, filename)
    end
end
function MISTBCGridv2(grid::AbstractString)
    if haskey(registry, grid)
        fname = mist_processed_fname(@datadep_str(grid))
    else
        grid_norm = normalize(grid; casefold=true)
        find_func = Base.Fix2(occursin, grid_norm)
        if find_func("niriss")
            fname = mist_processed_fname(datadep"MIST_NIRISS_v2.5")
        elseif find_func("jwst")
            fname = mist_processed_fname(datadep"MIST_JWST_v2.5")
        elseif find_func("wfpc2")
            fname = mist_processed_fname(datadep"MIST_HST_WFPC2_v2.5")
        elseif find_func("wfc3")
            fname = mist_processed_fname(datadep"MIST_HST_WFC3_v2.5")
        elseif mapreduce(find_func, &, ("acs", "hrc"))
            fname = mist_processed_fname(datadep"MIST_HST_ACS_HRC_v2.5")
        elseif mapreduce(find_func, &, ("acs", "sbc"))
            fname = mist_processed_fname(datadep"MIST_HST_ACS_SBC_v2.5")
        elseif mapreduce(find_func, &, ("acs", "wfc"))
            fname = mist_processed_fname(datadep"MIST_HST_ACS_WFC_v2.5")
        elseif find_func("hst")
            throw(ArgumentError("""Requested grid "$grid" unclear. Supported HST v2.5 grids are \
                                "acs_wfc", "acs_hrc", "acs_sbc", "wfpc2", and "wfc3"."""))
        elseif find_func("euclid")
            fname = mist_processed_fname(datadep"MIST_Euclid_v2.5")
        elseif find_func("lsst")
            fname = mist_processed_fname(datadep"MIST_LSST_v2.5")
        elseif find_func("wise")
            fname = mist_processed_fname(datadep"MIST_WISE_v2.5")
        elseif find_func("washington")
            fname = mist_processed_fname(datadep"MIST_Washington_v2.5")
        elseif find_func("swift")
            fname = mist_processed_fname(datadep"MIST_Swift_v2.5")
        elseif find_func("panstarrs")
            fname = mist_processed_fname(datadep"MIST_PanSTARRS_v2.5")
        elseif find_func("megacam")
            fname = mist_processed_fname(datadep"MIST_CFHT_MegaCam_v2.5")
        elseif find_func("skymapper")
            fname = mist_processed_fname(datadep"MIST_SkyMapper_v2.5")
        elseif find_func("uvit")
            fname = mist_processed_fname(datadep"MIST_UVIT_v2.5")
        elseif find_func("ukidss")
            fname = mist_processed_fname(datadep"MIST_UKIDSS_v2.5")
        elseif find_func("vista")
            fname = mist_processed_fname(datadep"MIST_VISTA_v2.5")
        elseif find_func("roboao")
            fname = mist_processed_fname(datadep"MIST_RoboAO_v2.5")
        elseif find_func("roman")
            fname = mist_processed_fname(datadep"MIST_Roman_v2.5")
        elseif mapreduce(find_func, |, ("johnson", "cousins", "bessell", "2mass", "kepler",
                                        "hipparcos", "tycho", "gaia", "tess"))
            fname = mist_processed_fname(datadep"MIST_UBVRIplus_v2.5")
        elseif find_func("hsc")
            fname = mist_processed_fname(datadep"MIST_HSC_v2.5")
        elseif find_func("spitzer")
            fname = mist_processed_fname(datadep"MIST_Spitzer_v2.5")
        elseif mapreduce(find_func, |, ("sdss", "sloan"))
            fname = mist_processed_fname(datadep"MIST_SDSS_v2.5")
        elseif find_func("galex")
            fname = mist_processed_fname(datadep"MIST_GALEX_v2.5")
        elseif find_func("iphas")
            fname = mist_processed_fname(datadep"MIST_IPHAS_v2.5")
        elseif mapreduce(find_func, |, ("splus", "s+"))
            fname = mist_processed_fname(datadep"MIST_SPLUS_v2.5")
        elseif find_func("decam")
            fname = mist_processed_fname(datadep"MIST_DECam_v2.5")
        else
            throw(ArgumentError("Unrecognized grid \"$grid\" for MIST v2.5."))
        end
    end
    return MISTBCGridv2(read_mist_bc_processed(fname), fname)
end
(grid::MISTBCGridv2)(feh::Real, afe::Real, Av::Real) = MISTBCTablev2(grid, feh, afe, Av)
Base.show(io::IO, z::MISTBCGridv2) =
    print(io, "MIST v2.5 bolometric correction grid for photometric system ",
          splitpath(splitext(z.filename)[1])[end])
Table(grid::MISTBCGridv2) = grid.table
Base.extrema(::Type{<:MISTBCGridv2}) =
    (lgTef = (first(gridinfov2.lgTef), last(gridinfov2.lgTef)),
     logg  = (first(gridinfov2.logg),  last(gridinfov2.logg)),
     feh   = (first(gridinfov2.feh),   last(gridinfov2.feh)),
     afe   = (first(gridinfov2.afe),   last(gridinfov2.afe)),
     Av    = (first(gridinfov2.Av),    last(gridinfov2.Av)),
     Rv    = (first(gridinfov2.Rv),    last(gridinfov2.Rv)))
filternames(grid::MISTBCGridv2) = columnnames(grid)[length(_mist_v2_dependents)+1:end]
gridname(::Type{<:MISTBCGridv2}) = "MIST"
zeropoints(::MISTBCGridv2) = zpt
chemistry(::Type{<:MISTBCGridv2}) = MISTChemistryv2()

######################################
# MISTBCTablev2
######################################

"""
    MISTBCTablev2

A MIST v2.5 bolometric correction table with fixed \\[Fe/H\\], \\[α/Fe\\], and ``A_V``.
Callable with `(Teff [K], logg [cgs])` to interpolate BCs.
"""
struct MISTBCTablev2{A <: Real, B, N} <: AbstractBCTable{A}
    feh::A
    afe::A
    Av::A
    itp::B
    filters::Tuple{Vararg{Symbol, N}}
end
MISTBCTablev2(feh::Real, afe::Real, Av::Real, itp, filters) =
    MISTBCTablev2(promote(feh, afe, Av)..., itp, filters)
Base.show(io::IO, z::MISTBCTablev2) =
    print(io, "MIST v2.5 bolometric correction table with [Fe/H] ", z.feh,
          ", [α/Fe] ", z.afe, ", and V-band extinction ", z.Av)
filternames(table::MISTBCTablev2) = table.filters
gridname(::Type{<:MISTBCTablev2}) = "MIST"
zeropoints(::MISTBCTablev2) = zpt
Base.extrema(::Type{<:MISTBCTablev2}) =
    (Teff = (10^first(gridinfov2.lgTef), 10^last(gridinfov2.lgTef)),
     logg = (first(gridinfov2.logg), last(gridinfov2.logg)))
MH(t::MISTBCTablev2) = t.feh
Z(t::MISTBCTablev2) = Z(chemistry(t), MH(t))
chemistry(::Type{<:MISTBCTablev2}) = MISTChemistryv2()

"""
    MISTBCTablev2(grid::MISTBCGridv2, feh::Real, afe::Real, Av::Real)

Interpolate the MIST v2.5 BCs in `grid` to fixed \\[Fe/H\\] (`feh`), \\[α/Fe\\] (`afe`),
and ``A_V`` (`Av`). Returns an instance callable with `(Teff [K], logg [cgs])`.

```jldoctest
julia> grid = MISTBCGridv2("GALEX")
MIST v2.5 bolometric correction grid for photometric system MIST_GALEX_v2.5

julia> table = MISTBCTablev2(grid, -1.0, 0.2, 0.0)
MIST v2.5 bolometric correction table with [Fe/H] -1.0, [α/Fe] 0.2, and V-band extinction 0.0

julia> length(table(5500, 4.5)) == 2
true
```
"""
function MISTBCTablev2(grid::MISTBCGridv2, feh::Real, afe::Real, Av::Real)
    check_vals_v2(feh, afe, Av)
    filters = filternames(grid)
    table = Table(grid)

    if feh ∈ _mist_v2_feh && afe ∈ _mist_v2_afe && Av ∈ _mist_v2_Av
        subtable = select_subtable_v2(table, feh, afe, Av)
        submatrix = Tables.matrix(getproperties(subtable, filters))
    else
        # Trilinear interpolation over (feh, afe, Av)
        feh_exact = feh ∈ _mist_v2_feh
        afe_exact = afe ∈ _mist_v2_afe
        Av_exact  = Av  ∈ _mist_v2_Av

        feh_idx = feh_exact ? findfirst(==(feh), _mist_v2_feh) :
                              searchsortedfirst(SVector(_mist_v2_feh), feh) - 1
        feh1 = _mist_v2_feh[feh_exact ? feh_idx : feh_idx]
        feh2 = feh_exact ? feh1 : _mist_v2_feh[feh_idx + 1]

        afe_idx = afe_exact ? findfirst(==(afe), _mist_v2_afe) :
                              searchsortedfirst(SVector(_mist_v2_afe), afe) - 1
        afe1 = _mist_v2_afe[afe_exact ? afe_idx : afe_idx]
        afe2 = afe_exact ? afe1 : _mist_v2_afe[afe_idx + 1]

        Av_idx_l = Av_exact ? findfirst(==(Av), _mist_v2_Av) :
                              searchsortedfirst(SVector(_mist_v2_Av), Av) - 1
        Av1 = _mist_v2_Av[Av_exact ? Av_idx_l : Av_idx_l]
        Av2 = Av_exact ? Av1 : _mist_v2_Av[Av_idx_l + 1]

        # Gather the (up to 8) corner matrices and reduce
        function _mat(f, a, v)
            Tables.matrix(getproperties(select_subtable_v2(table, f, a, v), filters))
        end
        if feh_exact && afe_exact
            submatrix = interp1d(Av, Av1, Av2, _mat(feh1, afe1, Av1), _mat(feh1, afe1, Av2))
        elseif feh_exact && Av_exact
            submatrix = interp1d(afe, afe1, afe2, _mat(feh1, afe1, Av1), _mat(feh1, afe2, Av1))
        elseif afe_exact && Av_exact
            submatrix = interp1d(feh, feh1, feh2, _mat(feh1, afe1, Av1), _mat(feh2, afe1, Av1))
        elseif feh_exact
            submatrix = interp2d(afe, Av, afe1, afe2, Av1, Av2,
                                 _mat(feh1, afe1, Av1), _mat(feh1, afe2, Av1),
                                 _mat(feh1, afe1, Av2), _mat(feh1, afe2, Av2))
        elseif afe_exact
            submatrix = interp2d(feh, Av, feh1, feh2, Av1, Av2,
                                 _mat(feh1, afe1, Av1), _mat(feh2, afe1, Av1),
                                 _mat(feh1, afe1, Av2), _mat(feh2, afe1, Av2))
        elseif Av_exact
            submatrix = interp2d(feh, afe, feh1, feh2, afe1, afe2,
                                 _mat(feh1, afe1, Av1), _mat(feh2, afe1, Av1),
                                 _mat(feh1, afe2, Av1), _mat(feh2, afe2, Av1))
        else
            # Full trilinear: interpolate feh first, then afe, then Av
            m_1_1_1 = _mat(feh1, afe1, Av1); m_2_1_1 = _mat(feh2, afe1, Av1)
            m_1_2_1 = _mat(feh1, afe2, Av1); m_2_2_1 = _mat(feh2, afe2, Av1)
            m_1_1_2 = _mat(feh1, afe1, Av2); m_2_1_2 = _mat(feh2, afe1, Av2)
            m_1_2_2 = _mat(feh1, afe2, Av2); m_2_2_2 = _mat(feh2, afe2, Av2)
            # Interpolate along feh
            m_f_1_1 = interp1d(feh, feh1, feh2, m_1_1_1, m_2_1_1)
            m_f_2_1 = interp1d(feh, feh1, feh2, m_1_2_1, m_2_2_1)
            m_f_1_2 = interp1d(feh, feh1, feh2, m_1_1_2, m_2_1_2)
            m_f_2_2 = interp1d(feh, feh1, feh2, m_1_2_2, m_2_2_2)
            # Interpolate along afe
            m_f_a_1 = interp1d(afe, afe1, afe2, m_f_1_1, m_f_2_1)
            m_f_a_2 = interp1d(afe, afe1, afe2, m_f_1_2, m_f_2_2)
            # Interpolate along Av
            submatrix = interp1d(Av, Av1, Av2, m_f_a_1, m_f_a_2)
        end
    end
    itp = interpolate((SVector(_mist_v2_lgTef), SVector(_mist_v2_logg)),
                      repack_submatrix(submatrix, length(_mist_v2_lgTef), length(_mist_v2_logg), filters),
                      Gridded(Linear()))
    itp = extrapolate(itp, Flat())
    return MISTBCTablev2(feh, afe, Av, itp, filters)
end
# Accept linear Teff in Kelvin, convert to log10 internally
(table::MISTBCTablev2)(Teff::Real, logg::Real) = table.itp(log10(Teff), logg)

###############################################################
# MIST v2.5 chemical mixture — Grevesse & Sauval (1998) / GS98
###############################################################

"""
    MISTChemistryv2()

Singleton struct representing the MIST v2.5 chemical mixture, based on [Grevesse1998](@cite) solar abundances (in contrast to [Asplund2009](@cite) used in v1.2).

Protosolar (initial) values [`X`](@ref), [`Y`](@ref), and [`Z`](@ref) are taken from the solar calibration in Table 1 of [Dotter2026](@cite). Photospheric values are derived from the
[GS98](@cite Grevesse1998) ratio `(Z/X)_phot = 0.023` and the calibrated surface helium abundance `Ysurf = 0.2482`. The primordial helium abundance [`Y_p`](@ref) is from [PlanckCollaboration2016](@cite).
"""
struct MISTChemistryv2 <: AbstractChemicalMixture end

X(::MISTChemistryv2)      = 0.7080          # Xinitial from Table 1, Dotter+2026
X_phot(::MISTChemistryv2) = 0.7518 / 1.023  # derived from GS98 (Z/X)=0.023, Ysurf=0.2482
Y(::MISTChemistryv2)      = 0.2735          # Yinitial
Y_phot(::MISTChemistryv2) = 0.2482          # Ysurf from Table 1, Basu 2004
Y_p(::MISTChemistryv2)    = 0.249           # Planck 2016
Z(::MISTChemistryv2)      = 0.0185          # Zinitial
Z_phot(::MISTChemistryv2) = 0.023 * (0.7518 / 1.023)  # ≈ 0.01690
function Y(mix::MISTChemistryv2, Zval)
    yp = Y_p(mix)
    return yp + ((Y(mix) - yp) / Z(mix)) * Zval
end
function MH(mix::MISTChemistryv2, Zval)
    Xval = X(mix, Zval)
    return log10((Zval / Xval) * (X(mix) / Z(mix)))
end
function Z(mix::MISTChemistryv2, MHval)
    solX, solY, solZ, yp = X(mix), Y(mix), Z(mix), Y_p(mix)
    zoverx = exp10(MHval + log10(solZ/solX))
    return zoverx * (1 - yp) / (1 + zoverx + zoverx * ((solY - yp) / solZ))
end