######################################
# MIST v1.2 grid constants
######################################

""" Unique values of `Teff` in the MIST v1.2 BC tables. """
const _mist_Teff = (2500.0, 2800.0, 3000.0, 3200.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0,
                    4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0, 6250.0, 6500.0, 6750.0,
                    7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0, 8500.0, 8750.0, 9000.0,
                    9250.0, 9500.0, 9750.0, 10000.0, 11000.0, 12000.0, 13000.0, 14000.0,
                    15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 25000.0, 30000.0,
                    35000.0, 40000.0, 45000.0, 50000.0, 60000.0, 70000.0, 80000.0, 90000.0,
                    100000.0, 110000.0, 120000.0, 130000.0, 140000.0, 150000.0, 160000.0,
                    170000.0, 180000.0, 190000.0, 200000.0, 300000.0, 400000.0, 500000.0,
                    600000.0, 700000.0, 800000.0, 900000.0, 1.0e6) # length=70
""" Unique values of `logg` in the MIST v1.2 BC tables. """
const _mist_logg = (-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
                    3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5) # length=26
""" Unique values of ``A_V`` in the MIST v1.2 BC tables. """
const _mist_Av = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0) # length=13
""" Unique values of ``R_V`` in the MIST v1.2 BC tables. """
const _mist_Rv = (3.1,) # length=1
""" Unique values of [Fe/H] in the MIST v1.2 BC tables. """
const _mist_feh = (-4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0,
                   -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75) # length=18

""" Unique values for dependent variables in the MIST v1.2 bolometric correction grid. """
const gridinfo = (Teff = _mist_Teff,
                  logg = _mist_logg,
                  feh = _mist_feh,
                  Av = _mist_Av,
                  Rv = _mist_Rv)
@compat public gridinfo

######################################
# MIST v1.2 validation and table I/O
######################################

"""
    check_vals(feh, Av)

Validate that [Fe/H] value `feh` and ``A_V`` value `Av` are valid for the MIST v1.2 BC grid.
Throws `ArgumentError` if check fails, returns `nothing` if check is successful.

```jldoctest
julia> using BolometricCorrections.MIST: check_vals

julia> check_vals(-2, 0.0) # Check passes, returns nothing

julia> using Test: @test_throws, Pass

julia> @test_throws(ArgumentError, check_vals(-5, 0.0)) isa Pass # Invalid `mh`, throws error
true

julia> @test_throws(ArgumentError, check_vals(-2, 100.0)) isa Pass # Invalid `Av`, throws error
true
```
"""
function check_vals(feh, Av)
    feh_ext = extrema(_mist_feh)
    if feh < first(feh_ext) || feh > last(feh_ext)
        throw(ArgumentError("Provided [Fe/H] $feh is outside the bounds for the MIST BC tables $feh_ext"))
    end
    Av_ext = extrema(_mist_Av)
    if Av < first(Av_ext) || Av > last(Av_ext)
        throw(ArgumentError("Provided A_v $Av is outside the bounds for the MIST BC tables $Av_ext"))
    end
end

# Extract a subtable out of table where table.feh == feh and table.Av == Av
function select_subtable(table::Table, feh::Real, Av::Real)
    @argcheck feh ∈ _mist_feh && Av ∈ _mist_Av
    # We can use the known structure of the data table to prevent having to do a filter at all
    # For each unique Rv (1), Av and feh there will be length(Teff) * length(logg) * length(Rv) rows
    nrows_per_Av  = length(_mist_Teff) * length(_mist_logg) * length(_mist_Rv) # 1820
    nrows_per_feh = nrows_per_Av * length(_mist_Av)
    idx1 = 1 + nrows_per_feh * (findfirst(==(feh), _mist_feh) - 1)
    idx2 = idx1 + nrows_per_Av * (findfirst(==(Av), _mist_Av) - 1)
    idx3 = idx2 + nrows_per_Av - 1
    return table[idx2:idx3]
end

######################################
# MISTBCGridv1
######################################

"""
    MISTBCGridv1(grid::AbstractString)

Load and return the MIST v1.2 bolometric corrections for the given photometric system `grid`.
This type is used to create instances of [`MISTBCTablev1`](@ref) with fixed dependent grid
variables (\\[Fe/H\\], Av). Call an instance with `(feh, Av)` or use the
[`MISTBCTablev1`](@ref) constructor directly.

```jldoctest
julia> grid = MISTBCGridv1("JWST")
MIST v1.2 bolometric correction grid for photometric system MIST_JWST

julia> grid(-1.01, 0.11) # Can be called to construct table with interpolated [Fe/H], Av
MIST v1.2 bolometric correction table with [Fe/H] -1.01 and V-band extinction 0.11
```
"""
struct MISTBCGridv1{A, B, C <: AbstractString} <: AbstractBCGrid{A}
    table::B # Usually a TypedTables.Table
    filename::C
    function MISTBCGridv1(table::B, filename::C) where {B, C}
        A = Base.promote_eltype(first(table))
        new{A, B, C}(table, filename)
    end
end
function MISTBCGridv1(grid::AbstractString)
    if haskey(registry, grid)
        fname = mist_processed_fname(@datadep_str(grid))
    else
        grid = normalize(grid; casefold=true)
        find_func = Base.Fix2(occursin, grid)
        if find_func("jwst")
            fname = mist_processed_fname(datadep"MIST_JWST")
        elseif find_func("wfpc2")
            fname = mist_processed_fname(datadep"MIST_HST_WFPC2")
        elseif find_func("wfc3")
            fname = mist_processed_fname(datadep"MIST_HST_WFC3")
        elseif mapreduce(find_func, &, ("acs", "hrc"))
            fname = mist_processed_fname(datadep"MIST_HST_ACS_HRC")
        elseif mapreduce(find_func, &, ("acs", "wfc"))
            fname = mist_processed_fname(datadep"MIST_HST_ACS_WFC")
        elseif find_func("hst")
            throw(ArgumentError("""Requested grid "$grid" unclear. Supported HST grids are \
                                "acs_wfc", "acs_hrc", "wfpc2", and "wfc3"."""))
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
        elseif mapreduce(find_func, |, ("johnson", "cousins", "bessell", "2mass", "kepler",
                                        "hipparcos", "tycho", "gaia", "tess"))
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
        elseif mapreduce(find_func, |, ("wfirst", "roman"))
            @info """The WFIRST (now the Nancy Grace Roman Telescope) bolometric corrections are based \
                     on a preliminary filter set from 2018."""
            fname = mist_processed_fname(datadep"MIST_WFIRST")
        elseif find_func("decam")
            fname = mist_processed_fname(datadep"MIST_DECam")
        else
            throw(ArgumentError("Unrecognized grid \"$grid\" for MIST v1.2."))
        end
    end
    return MISTBCGridv1(read_mist_bc_processed(fname), fname)
end
(grid::MISTBCGridv1)(feh::Real, Av::Real) = MISTBCTablev1(grid, feh, Av)
Base.show(io::IO, z::MISTBCGridv1) =
    print(io, "MIST v1.2 bolometric correction grid for photometric system ",
          splitpath(splitext(z.filename)[1])[end])
Table(grid::MISTBCGridv1) = grid.table
Base.extrema(::Type{<:MISTBCGridv1}) = (Teff = (first(gridinfo.Teff), last(gridinfo.Teff)),
                                        logg = (first(gridinfo.logg), last(gridinfo.logg)),
                                        feh  = (first(gridinfo.feh),  last(gridinfo.feh)),
                                        Av   = (first(gridinfo.Av),   last(gridinfo.Av)),
                                        Rv   = (first(gridinfo.Rv),   last(gridinfo.Rv)))
filternames(grid::MISTBCGridv1) = columnnames(grid)[length(_mist_v1_dependents)+1:end]
gridname(::Type{<:MISTBCGridv1}) = "MIST"
zeropoints(::MISTBCGridv1) = zpt
chemistry(::Type{<:MISTBCGridv1}) = MISTChemistryv1()

######################################
# MISTBCTablev1
######################################

"""
    MISTBCTablev1

A MIST v1.2 bolometric correction table with fixed \\[Fe/H\\] and V-band extinction `Av`.
Callable with arguments `(Teff [K], logg [cgs])` to interpolate BCs.
"""
struct MISTBCTablev1{A <: Real, B, N} <: AbstractBCTable{A}
    feh::A
    Av::A
    itp::B
    filters::Tuple{Vararg{Symbol, N}}
end
MISTBCTablev1(feh::Real, Av::Real, itp, filters) =
    MISTBCTablev1(promote(feh, Av)..., itp, filters)
Base.show(io::IO, z::MISTBCTablev1) =
    print(io, "MIST v1.2 bolometric correction table with [Fe/H] ", z.feh,
          " and V-band extinction ", z.Av)
filternames(table::MISTBCTablev1) = table.filters
gridname(::Type{<:MISTBCTablev1}) = "MIST"
zeropoints(::MISTBCTablev1) = zpt
Base.extrema(::Type{<:MISTBCTablev1}) =
    (Teff = (first(gridinfo.Teff), last(gridinfo.Teff)),
     logg = (first(gridinfo.logg), last(gridinfo.logg)))
MH(t::MISTBCTablev1) = t.feh
Z(t::MISTBCTablev1) = Z(chemistry(t), MH(t))
chemistry(::Type{<:MISTBCTablev1}) = MISTChemistryv1()

"""
    MISTBCTablev1(grid::MISTBCGridv1, feh::Real, Av::Real)

Interpolate the MIST v1.2 BCs in `grid` to fixed \\[Fe/H\\] (`feh`) and ``A_V`` (`Av`),
returning an instance callable with `(Teff [K], logg [cgs])`.

```jldoctest
julia> grid = MISTBCGridv1("JWST")
MIST v1.2 bolometric correction grid for photometric system MIST_JWST

julia> table = MISTBCTablev1(grid, -1.01, 0.011)
MIST v1.2 bolometric correction table with [Fe/H] -1.01 and V-band extinction 0.011

julia> length(table(2755, 0.01)) == 29
true

julia> size(table([2755, 2756], [0.01, 0.02]))
(29, 2)

julia> using TypedTables: Table

julia> table(Table, [2755, 2756], [0.01, 0.02]) isa Table
true
```
"""
function MISTBCTablev1(grid::MISTBCGridv1, feh::Real, Av::Real)
    check_vals(feh, Av)
    filters = filternames(grid)
    table = Table(grid)

    if feh ∈ _mist_feh && Av ∈ _mist_Av
        subtable = select_subtable(table, feh, Av)
        submatrix = Tables.matrix(getproperties(subtable, filters))
    else
        if feh ∈ _mist_feh
            Av_idx = searchsortedfirst(SVector(_mist_Av), Av) - 1
            Av1, Av2 = _mist_Av[Av_idx], _mist_Av[Av_idx + 1]
            mat1 = Tables.matrix(getproperties(select_subtable(table, feh, Av1), filters))
            mat2 = Tables.matrix(getproperties(select_subtable(table, feh, Av2), filters))
            submatrix = interp1d(Av, Av1, Av2, mat1, mat2)
        elseif Av ∈ _mist_Av
            feh_idx = searchsortedfirst(SVector(_mist_feh), feh) - 1
            feh1, feh2 = _mist_feh[feh_idx], _mist_feh[feh_idx + 1]
            mat1 = Tables.matrix(getproperties(select_subtable(table, feh1, Av), filters))
            mat2 = Tables.matrix(getproperties(select_subtable(table, feh2, Av), filters))
            submatrix = interp1d(feh, feh1, feh2, mat1, mat2)
        else
            feh_idx = searchsortedfirst(SVector(_mist_feh), feh) - 1
            feh1, feh2 = _mist_feh[feh_idx], _mist_feh[feh_idx + 1]
            Av_idx = searchsortedfirst(SVector(_mist_Av), Av) - 1
            Av1, Av2 = _mist_Av[Av_idx], _mist_Av[Av_idx + 1]
            mat1_1 = Tables.matrix(getproperties(select_subtable(table, feh1, Av1), filters))
            mat2_1 = Tables.matrix(getproperties(select_subtable(table, feh2, Av1), filters))
            mat1_2 = Tables.matrix(getproperties(select_subtable(table, feh1, Av2), filters))
            mat2_2 = Tables.matrix(getproperties(select_subtable(table, feh2, Av2), filters))
            submatrix = interp2d(feh, Av, feh1, feh2, Av1, Av2, mat1_1, mat2_1, mat1_2, mat2_2)
        end
    end
    itp = interpolate((SVector(_mist_logg), SVector(_mist_Teff)),
                      repack_submatrix(submatrix, length(_mist_logg), length(_mist_Teff), filters),
                      Gridded(Linear()))
    itp = extrapolate(itp, Flat())
    return MISTBCTablev1(feh, Av, itp, filters)
end
(table::MISTBCTablev1)(Teff::Real, logg::Real) = table.itp(logg, Teff)

##############################
# Chemical mixture information

"""
    MISTChemistryv1()
Returns a singleton struct representing the MIST v1.2 chemical mixture model.
MIST v1.2 assumes the *protostellar* [Asplund2009](@citet) solar abundances. Sum of protostellar hydrogen, helium,
metal mass fractions from last row of Table 4 sums to 0.9999, not 1 as it should.
To keep calculations consistent, the protostellar values are normalized to sum to 1 here.

```jldoctest
julia> using BolometricCorrections.MIST: MISTChemistryv1, X, Y, Z, X_phot, Y_phot, Z_phot,
                                         MH;

julia> chem = MISTChemistryv1();

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
struct MISTChemistryv1 <: AbstractChemicalMixture end
X(::MISTChemistryv1) = 0.7154 / 0.9999
X_phot(::MISTChemistryv1) = 0.7381
Y(::MISTChemistryv1) = 0.2703 / 0.9999
Y_phot(::MISTChemistryv1) = 0.2485
Y_p(::MISTChemistryv1) = 0.249
Z(::MISTChemistryv1) = 0.0142 / 0.9999
Z_phot(::MISTChemistryv1) = 0.0134
function Y(mix::MISTChemistryv1, Zval)
    yp = Y_p(mix)
    return yp + ((Y(mix) - yp) / Z(mix)) * Zval
end
function MH(mix::MISTChemistryv1, Zval)
    Xval = X(mix, Zval)
    solZ = Z(mix)
    solX = X(mix)
    # return log10(Zval / Xval) - log10(solZ / solX)
    # Fuse expression into one log10 call for efficiency
    return log10((Zval / Xval) * (solX / solZ))
end
function Z(mix::MISTChemistryv1, MHval)
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
