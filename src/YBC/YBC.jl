module YBC

# using ..BolometricCorrections: @compat
using ..BolometricCorrections: AbstractChemicalMixture, AbstractBCGrid, AbstractBCTable, AbstractMassLoss, Bjorklund2021MassLoss, AllHardwareNumeric, without, interp2d, interp1d
import ..BolometricCorrections: X, X_phot, Y, Y_phot, Y_p, Z, Z_phot, MH, chemistry, filternames, zeropoints

using Compat: @compat
import CSV
using StaticArrays: SVector, sacollect
using TypedTables: Table

"""YBC data are stored in FITS files with natural datatype Float32."""
const dtype = Float32
const HardwareNumeric = without(dtype, AllHardwareNumeric)

# Code to initialize data storage mechanisms
include("init.jl")



"""
    parse_filterinfo(f::AbstractString)
Given a `filter.info` file from YBC, return a `TypedTables.Table` with the corresponding information. The columns are described as follows by YBC:

`filter.info`: infomation for the filters. The columns are filter number, filter name, filter file used, effective wavelength (in Angstrom), band width, absolute flux of Vega (in unit of erg cm-2 s-1 A-1), photometric system (AB, Vega, ST or TG), type of device (energy counting or number counting for CCD), reference magnitude of Vega applied, irrelevant value, comments.
"""
function parse_filterinfo(f::AbstractString)
    # Contains final column with # <comments>, we don't want these so filter first
    cleaned = IOBuffer(join(map(line -> split(line, "#")[1], readlines(String(f))), "\n"))
    return CSV.read(cleaned, Table; comment="#", delim=' ', ignorerepeated=true, 
        header=["index", "names", "file", "effective_wavelength", "width", "flux_zeropoint", "photometric_system", "detector_type", "mag_zeropoint", "blank"])
    # types=[Int, String, String, dtype, dtype, dtype, String, String, dtype, Int, String])
end

# Presently only supporting the standard YBC BCs, non-rotating
check_prefix(prefix) = prefix != "YBC" ? throw(ArgumentError("""prefix = $prefix not presently supported -- use "YBC".""")) : return nothing

"""
    check_vals(mh, Av, gridinfo::NamedTuple)

Validate that [M/H] value `mh` and ``A_V`` value `Av` are valid for the YBC with grid information `gridinfo`.
This function expects the available [M/H] values to be `gridinfo.MH` and ``A_V`` values to be `gridinfo.Av`.
Throws `ArgumentError` if check fails, returns `nothing` if check is successful.

```jldoctest
julia> using BolometricCorrections.YBC.PHOENIX: check_vals, gridinfo

julia> check_vals(-2, 0.0, gridinfo) # Check passes, returns nothing

julia> using Test: @test_throws, Pass

julia> @test_throws(ArgumentError, check_vals(-5, 0.0, gridinfo)) isa Pass # Invalid `mh`, throws error
true

julia> @test_throws(ArgumentError, check_vals(-2, 100.0, gridinfo)) isa Pass # Invalid `Av`, throws error
true
```
"""
function check_vals(mh, Av, gridinfo::NamedTuple)
    mh_ext = extrema(gridinfo.MH)
    if mh < first(mh_ext) || mh > last(mh_ext)
        throw(ArgumentError("Provided [M/H] $mh is outside the bounds for the BC grid $mh_ext"))
    end
    Av_ext = extrema(gridinfo.Av)
    if Av < first(Av_ext) || Av > last(Av_ext)
        throw(ArgumentError("Provided A_v $Av is outside the bounds for the BC grid $Av_ext"))
    end
end

##########################################################################

"""
    PARSECChemistry()
Returns a singleton struct representing the PARSEC chemical mixture model.
We presently only include scaled-solar models. The solar protostellar chemical
mixture for PARSEC was calibrated to reproduce solar photospheric observations
via a forward modeling approach (see section 3 of [Bressan2012](@citet)). The
full solar calibration assumed for PARSEC is given in Table 3 of [Bressan2012](@citet).
The distribution of heavy metals is taken from [Grevesse1998](@citet) and [Caffau2011](@citet) (see section 4 of [Bressan2012](@citet)).

```jldoctest
julia> using BolometricCorrections.YBC: PARSECChemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

julia> chem = PARSECChemistry();

julia> X(chem) + Y(chem) + Z(chem) ≈ 1 # solar protostellar values
true

julia> X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1 # solar photospheric values
true

julia> MH(chem, Z(chem) * 0.1) ≈ -0.9400696788068212
true

julia> Z(chem, -0.9400696788068212) ≈ Z(chem) * 0.1
true
```
"""
struct PARSECChemistry <: AbstractChemicalMixture end
# For PARSEC, a choice can be made as to whether the initial solar
# chemical composition is taken to be the observed reference value
# i.e., Z⊙, Y⊙ in Table 3 of Bressan2012, or the photospheric abundances
# of the best-fit solar calibration model, i.e., Zs, Ys in Table 3. 
# For consistency with PARSEC's conversion between Z and [M/H], we will
# assume the observed reference values.

X(mix::PARSECChemistry) = 1 - Y(mix) - Z(mix) # 0.70226
X_phot(mix::PARSECChemistry) = 1 - Y_phot(mix) - Z_phot(mix)  # 0.73616
Y(::PARSECChemistry) = 0.28 # Y_initial in Table 3 of Bressan2012
Y_phot(::PARSECChemistry) = 0.2485  # Y⊙ in Table 3 of Bressan2012
# Y_phot(::PARSECChemistry) = 0.24787 # Y_S in Table 3 of Bressan2012
Z(::PARSECChemistry) = 0.01774 # Z_initial in Table 3 of Bressan2012
Z_phot(::PARSECChemistry) = 0.01524 # 0.01774 # Z⊙ in Table 3 of Bressan2012
# Z_phot(::PARSECChemistry) = 0.01597 # Z_S in Table 3 of Bressan2012
Y_p(::PARSECChemistry) = 0.2485

Y(mix::PARSECChemistry, Zval) = Y_p(mix) + 178//100 * Zval # γ = 1.78
# X generic
MH(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval)) - log10(Z_phot(mix) / X_phot(mix))
# MH(mix::PARSECChemistry, Zval) = log10(Zval / X(mix, Zval) / Z(mix) * X(mix))
function Z(mix::PARSECChemistry, MHval)
    # [M/H] = log(Z/X) - log(Z/X)☉ with Z☉ = solz
    # Z/X = exp10( [M/H] + log(Z/X)☉ )
    # X = 1 - Y - Z
    # Y ≈ Y_p + γ * Z for parsec (see Y(mix::PARSECChemistry, Zval) above)
    # so X ≈ 1 - (Y_p + γ * Z) - Z = 1 - Y_p - (1 + γ) * Z
    # Substitute into line 2,
    # Z / (1 - Y_p - (1 + γ) * Z) = exp10( [M/H] + log(Z/X)☉ )
    # Z = (1 - Y_p - (1 + γ) * Z) * exp10( [M/H] + log(Z/X)☉ )
    # let A = exp10( [M/H] + log(Z/X)☉ )
    # Z = (1 - Y_p) * A - (1 + γ) * Z * A
    # Z + (1 + γ) * Z * A = (1 - Y_p) * A
    # Z (1 + (1 + γ) * A) = (1 - Y_p) * A
    # Z = (1 - Y_p) * A / (1 + (1 + γ) * A)
    zoverx = exp10(MHval + log10(Z_phot(mix) / X_phot(mix)))
    γ = 178//100
    return (1 - Y_p(mix)) * zoverx / (1 + (1 + γ) * zoverx)
end

##########################################################################

include("PHOENIX.jl")
using .PHOENIX
@compat public PHOENIX
export PHOENIXYBCTable, PHOENIXYBCGrid

include("ATLAS9.jl")
using .ATLAS9
@compat public ATLAS9
export ATLAS9YBCTable, ATLAS9YBCGrid

include("Koester.jl")
using .KoesterWD
# Koester should not be public, don't export, only used internally

include("WMbasic.jl")
using .WMbasic
# WMBasic should not be public, don't export, only used internally

##########################################################################

"""
    YBCGrid(grid::AbstractString;
            mass_loss_model::AbstractMassLoss = Bjorklund2021MassLoss(),
            phoenix_transition=(Teff=(5500f0, 6500f0), logg=(-Inf32, Inf32)), # Lower than Teff[1], use PHOENIX
            wmbasic_transition=(Teff=(18000f0, 20000f0), logg=(2.5f0, 3.5f0), Mdot=(-Inf32, 1f-8)),
            koester_transition=(Teff=(5800f0, 6300f0), logg=(5.0f0, 6.0f0)))

Load and return the [YBC](@citet Chen2019) bolometric corrections for the given photometric system `grid`,
which must be a valid entry in `BolometricCorrections.YBC.systems`.
This model interpolates between bolometric correction grids derived from
several different atmoshere libraries -- the individual grids that make up
the integrated `YBCGrid` are [PHOENIX](@ref YBCPHOENIX), [ATLAS9](@ref YBCATLAS9),
the [Koester & Tremblay white dwarf grid](@ref YBCKoester), and the 
[WM-basic O and B star grid](@ref YBCWMbasic).

This type is used to create instances of [`YBCTable`](@ref) that have fixed dependent
grid variables (\\[M/H\\], Av). This can be done either by calling an instance of
`YBCGrid` with `(mh, Av)` arguments or by using the appropriate constructor for [`YBCTable`](@ref).

```jldoctest
julia> using BolometricCorrections.YBC: YBCGrid

julia> grid = YBCGrid("acs_wfc")
YBC WM-basic bolometric correction grid for photometric system YBC/acs_wfc.

julia> grid(-1.01, 0.11) # Can be called to construct table with interpolated [M/H], Av
YBC bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.11
```
"""
# struct YBCGrid{AA <: Number, A1, A2, A3, A4, B, C <: AbstractMassLoss, D <: AbstractVector{AA}, N} <: AbstractBCGrid{AA}
struct YBCGrid{AA <: Number, A, B, C <: AbstractMassLoss, D <: AbstractVector{AA}, N} <: AbstractBCGrid{AA}
    grids::A
    # phoenix::A1
    # atlas9::A2
    # wmbasic::A3
    # koester::A4
    limits::B
    mass_loss_model::C
    mag_zpt::D
    systems::Vector{String}
    name::String
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end

function YBCGrid(grids, limits, mass_loss_model, mag_zpt, systems, name, filternames)
    return YBCGrid(grids, limits, mass_loss_model, mag_zpt, String.(systems), String(name), tuple(Symbol.(filternames)...))
end
# function YBCGrid(phoenix, atlas9, wmbasic, koester, limits, mass_loss_model, mag_zpt, systems, name, filternames)
#     return YBCGrid(phoenix, atlas9, wmbasic, koester, limits, mass_loss_model, mag_zpt, String.(systems), String(name), tuple(Symbol.(filternames)...))
# end

function YBCGrid(grid::AbstractString;
                 mass_loss_model::AbstractMassLoss = Bjorklund2021MassLoss(),
                 prefix::AbstractString="YBC",
                #  phoenix_transition=(Teff=(5500f0, 6500f0),), # logg=(-Inf32, Inf32)), # Lower than Teff[1], use PHOENIX
                #  wmbasic_transition=(Teff=(18000f0, 20000f0), logg=(2.5f0, 3.5f0), Mdot=(-Inf32, 1f-8)),
                #  koester_transition=(Teff=(5800f0, 6300f0), logg=(5.0f0, 6.0f0)))
                 phoenix_lim=(Teff=(0f0, 6500f0), logg=(-Inf32, 6.0f0), Mdot=(0f0, 2f-8)),
                 atlas9_lim=(Teff=(5500f0, Inf32), logg=(-Inf32, 5.5f0), Mdot=(0f0, 2f-8)),
                 wmbasic_lim=(Teff=(18000f0, Inf32), logg=(-Inf32, 5.5f0), Mdot=(1f-8, Inf32)),
                 koester_lim=(Teff=(5000f0, Inf32), logg=(5.75f0, Inf32), Mdot=(0f0, 2f-8)))

    check_prefix(prefix)
    # Grid order is assumed later; this *cannot* be safely changed on its own
    grids = (phoenix = PHOENIXYBCGrid(grid; prefix=prefix),
             atlas9 = ATLAS9YBCGrid(grid; prefix=prefix),
             wmbasic = WMbasic.WMbasicYBCGrid(grid; prefix=prefix),
             koester = KoesterWD.KoesterWDYBCGrid(grid; prefix=prefix))
    limits = (phoenix = phoenix_lim, atlas9 = atlas9_lim, wmbasic = wmbasic_lim, koester = koester_lim)
    f = first(grids)
    return YBCGrid(grids, limits, mass_loss_model, f.mag_zpt, f.systems, f.name, f.filters)
end
(grid::YBCGrid)(mh::Real, Av::Real) = YBCTable(grid, mh, Av)
Base.show(io::IO, z::YBCGrid) = print(io, "YBC bolometric correction grid for photometric system $(z.name).")
filternames(grid::YBCGrid) = grid.filters
chemistry(::Type{<:YBCGrid}) = PARSECChemistry()
const _ybc_extrema = begin
    function _map_extrema(key::Symbol, ex)
        return (minimum(getproperty(i, key)[1] for i in ex if key in keys(i)), maximum(getproperty(i, key)[2] for i in ex if key in keys(i)))
    end
    ex = [extrema(g) for g in (PHOENIXYBCGrid, ATLAS9YBCGrid, KoesterWD.KoesterWDYBCGrid, WMbasic.WMbasicYBCGrid)]
    kk = (:Teff, :logg, :MH, :Av, :Mdot, :Rv)
    NamedTuple{kk}(_map_extrema(k, ex) for k in kk)
end
Base.extrema(::Type{<:YBCGrid}) = _ybc_extrema
# function Table(grid::WMbasicYBCGrid)
#     data = grid.data
#     tables = Vector{Table}(undef, length(data))
# end

"""
    YBCTable(grid::YBCGrid, mh::Real, Av::Real)

Interpolates the YBC bolometric corrections in `grid` to a fixed value of \\[M/H\\]
(`mh`), V-band extinction (`Av`), leaving only `Teff`, `logg`, and `Mdot` as dependent
variables. Returns an instance that is callable with arguments `(Teff [K], logg [cgs], 
Mdot [solMass / yr])` to interpolate the bolometric corrections as a function
of temperature, surface gravity, and mass outflow rate.

```jldoctest ybctable
julia> using BolometricCorrections.YBC: YBCGrid, YBCTable

julia> grid = YBCGrid("acs_wfc")
YBC bolometric correction grid for photometric system YBC/acs_wfc.

julia> table = YBCTable(grid, -1.01, 0.011) # Interpolate table from full grid
YBC bolometric correction table for system YBC/acs_wfc with [M/H] -1.01 and V-band extinction 0.011

julia> length(table(25_0254.0, 2.54, 5e-6)) == 12 # Returns BC in each filter
true

julia> size(table([25_0254.0, 25_0354.0], [2.54, 2.56], [4e-6, 5e-6])) # `table(array, array)` is also supported
(12, 2)

julia> using TypedTables: Table # `table(Table, array, array)` will return result as a Table

julia> table(Table, [25_0254.0, 25_0354.0], [2.54, 2.56], [4e-6, 5e-6]) isa Table
true
```

The YBC [WM-basic](@ref YBCWMbasic) BCs require the mass outflow rate `Mdot` -- when called with two arguments, 
it is assumed the arguments are `Teff, logg` and that `Mdot=0` such that the WM-basic models
are not used. 

```jldoctest ybctable
julia> table(25_0254.0, 2.54) == table(25_0254.0, 2.54, 0.0)
true
```

If called with via the one-argument method, taking types like `NamedTuple`,
required parameters will be parsed from the provided argument. If quantities sufficient
to estimate a mass-loss rate are identified, then the mass-loss rate will be automatically
estimated and the WM-basic models will be used when appropriate. An example of this usage 
is given below -- the default [`Bjorklund2021MassLoss`](@ref) model is used,
which calculates the stellar mass-loss rate from the metallicity (i.e., `Z(table)`)
and the luminosity `logL`. `logg` and `Teff` are treated normally.

```jldoctest ybctable
julia> table((logg = 2.54, Teff = 25_054.0, logL = 5)) isa AbstractVector
true
```
"""
struct YBCTable{A <: Real, B, C, D <: AbstractMassLoss, N} <: AbstractBCTable{A}
    MH::A
    Av::A
    mag_zpt::Vector{A}
    systems::Vector{String}
    name::String
    tables::B
    limits::C
    mass_loss_model::D
    filters::Tuple{Vararg{Symbol, N}} # NTuple{N, Symbol} giving filter names
end
function YBCTable(MH::Real, Av::Real, mag_zpt::Vector{<:Real}, systems, name, tables, limits, mass_loss_model::AbstractMassLoss, filters)
    T = dtype # promote_type(typeof(MH), typeof(Av), eltype(mag_zpt))
    return YBCTable(convert(T, MH), convert(T, Av), convert(Vector{T}, mag_zpt), convert.(String, systems), String(name), tables, limits, mass_loss_model, filters)
end
Base.show(io::IO, z::YBCTable) = print(io, "YBC bolometric correction table for system $(z.name) with [M/H] ",
                                       MH(z), " and V-band extinction ", z.Av)
filternames(table::YBCTable) = table.filters
# zeropoints(table::YBCTable) = table.mag_zpt
chemistry(::Type{<:YBCTable}) = PARSECChemistry()
MH(t::YBCTable) = t.MH
Z(t::YBCTable) = Z(chemistry(t), MH(t))
function Base.extrema(::Type{<:YBCTable})
    ex = _ybc_extrema
    return (Teff = ex.Teff, logg = ex.logg, Mdot = ex.Mdot)
end

function YBCTable(grid::YBCGrid, mh::Real, Av::Real)
    check_vals(mh, Av, extrema(grid))
    filters = filternames(grid)

    # tables = [g(mh, Av) for g in grid.grids]
    tables = NamedTuple{keys(grid.grids)}(g(mh, Av) for g in grid.grids)
    return YBCTable(mh, Av, grid.mag_zpt, grid.systems, grid.name, tables, grid.limits, grid.mass_loss_model, filters)
end

# Evaluation methods
(table::YBCTable)(Teff::Real, logg::Real) = table(Teff, logg, zero(dtype))
"""
    check_vals(nt::NamedTuple, Teff, logg, Mdot)
Check if `(Teff, logg, Mdot)` all fall within the bounds set in `nt`, which might look like `nt = (Teff=(0f0, 6500f0), logg=(-Inf32, 6.0f0), Mdot=(0f0, 2f-8))`."""
check_vals(nt::NamedTuple, Teff, logg, Mdot) = (first(nt.Teff) <= Teff <= last(nt.Teff)) && (first(nt.logg) <= logg <= last(nt.logg)) && (first(nt.Mdot) <= Mdot <= last(nt.Mdot))
"""
    _get_table(tables, i::Int)
Improve type stability in `(table::YBCTable)(Teff::Real, logg::Real, Mdot::Real)` by converting indexing calls like `table.tables[1]` into the corresponding NamedTuple keys like `tables.phoenix`. Assumes a particular mapping between index and key: `((i=1, key=:phoenix), (i=2, key=:atlas9), (i=3, key=:wmbasic), (i=4, key=:koester))`. This function is type unstable but not terribly slow."""
function _get_table(tables, i::Int)
# function _get_table(tables, limits, i::Int)
    if i == 1
        # return tables.phoenix, limits.phoenix # ::PHOENIXYBCTable
        return tables.phoenix
        # return :phoenix
    elseif i == 2
        # return tables.atlas9, limits.atlas9 # ::ATLAS9YBCTable
        return tables.atlas9
        # return :atlas9
    elseif i == 3
        # return tables.wmbasic, limits.wmbasic # ::WMbasic.WMbasicYBCTable
        return tables.wmbasic
        # return :wmbasic
    elseif i == 4
        # return tables.koester, limits.koester # ::KoesterWD.KoesterWDYBCTable
        return tables.koester
        # return :koester
    else
        error("Invalid index $i; must be between 1 and 4.")
    end
end
# function _get_limits(tables, limits, i::Int)
#     if i == 1
#         return tables.phoenix, limits.phoenix # ::PHOENIXYBCTable
#         # return :phoenix
#     elseif i == 2
#         return tables.atlas9, limits.atlas9 # ::ATLAS9YBCTable
#         # return :atlas9
#     elseif i == 3
#         return tables.wmbasic, limits.wmbasic # ::WMbasic.WMbasicYBCTable
#         # return :wmbasic
#     elseif i == 4
#         return tables.koester, limits.koester # ::KoesterWD.KoesterWDYBCTable
#         # return :koester
#     else
#         error("Invalid index $i; must be between 1 and 4.")
#     end
# end
"""
    _itp_tables(table1, table2, lim1, lim2, Teff, logg, Mdot)
Given two YBC subtables `(table1, table)` that are to be used over the limits `(lim1, lim2)`, interpolate between them at `(Teff, logg, Mdot)` with weighting based on how close the point is to the boundaries defined in `lim`.
"""
function _itp_tables(table1, table2, lim1, lim2, Teff, logg, Mdot)
    # Evaluate table1, table2
    arg = (Teff = Teff, logg = logg, Mdot = Mdot)
    r1, r2 = table1(arg), table2(arg)

    # Determine interpolation weights
    d_t1 = (abs(Teff - lim1.Teff[1]), abs(Teff - lim1.Teff[2])) # ul_t1 = ifelse(abs(Teff - lim1.Teff[1]) < abs(Teff - lim1.Teff[2]), 1, 2)
    ul_t1 = ifelse(d_t1[1] < d_t1[2], 1, 2)
    d_t2 = (abs(Teff - lim2.Teff[1]), abs(Teff - lim2.Teff[2]))
    ul_t2 = ifelse(d_t2[1] < d_t2[2], 1, 2) # ul_t2 = ifelse(abs(Teff - lim2.Teff[1]) < abs(Teff - lim2.Teff[2]), 1, 2)
    # logTe = log10(Teff)
    # d_t1 = (abs(logTe - log10(lim1.Teff[1])), abs(logTe - log10(lim1.Teff[2]))) # ul_t1 = ifelse(abs(Teff - lim1.Teff[1]) < abs(Teff - lim1.Teff[2]), 1, 2)
    # ul_t1 = ifelse(d_t1[1] < d_t1[2], 1, 2)
    # d_t2 = (abs(logTe - log10(lim2.Teff[1])), abs(logTe - log10(lim2.Teff[2])))
    # ul_t2 = ifelse(d_t2[1] < d_t2[2], 1, 2) # ul_t2 = ifelse(abs(Teff - lim2.Teff[1]) < abs(Teff - lim2.Teff[2]), 1, 2)

    d_g1 = (abs(logg - lim1.logg[1]), abs(logg - lim1.logg[2]))
    ul_g1 = ifelse(d_g1[1] < d_g1[2], 1, 2)
    d_g2 = (abs(logg - lim2.logg[1]), abs(logg - lim2.logg[2]))
    ul_g2 = ifelse(d_g2[1] < d_g2[2], 1, 2)
    # We'll determine the relative weights based on the log(Te) and log(g) differentials
    # Deterine if the limits for the tables are the same in any dimensions
    Teff_lim, logg_lim = (lim1.Teff[ul_t1], lim2.Teff[ul_t2]), (lim1.logg[ul_g1], lim2.logg[ul_g2])
    # t_eq, g_eq = (lim1.Teff[ul_t1] == lim2.Teff[ul_t2], lim1.logg[ul_g1] == lim2.logg[ul_g2])
    t_eq, g_eq = (Teff_lim[1] == Teff_lim[2], logg_lim[1] == logg_lim[2])
    if t_eq & g_eq
        # If both bounds equal, just return mean of two results
        return @. (r1 + r2) / 2
    elseif t_eq # lim1.Teff[ul_t1] == lim2.Teff[ul_t2]
        println("t_eq")
        # If temperature bounds equal, only interpolate in logg
        return interp1d(logg, lim1.logg[ul_g1], lim2.logg[ul_g2], r1, r2)
    elseif g_eq # lim1.logg[ul_g1] == lim2.logg[ul_g2]
        # If logg bounds equal, only interpolate in log(Teff)
        println("g_eq")
        if all(<(Teff), Teff_lim) # All Teff_lim < Teff, so we are not "between" two tables -- just return greater table
            return ifelse(Teff_lim[1] > Teff_lim[2], r1, r2)
        end
        # if all(<(logg), logg_lim)
        return interp1d(log10(Teff), log10(lim1.Teff[ul_t1]), log10(lim2.Teff[ul_t2]), r1, r2)
        # return interp1d(Teff, lim1.Teff[ul_t1], lim2.Teff[ul_t2], r1, r2)
    else
        # If no bounds equal, use 2D interpolation
        return interp2d(log10(Teff), logg, log10(lim1.Teff[ul_t1]), log10(lim2.Teff[ul_t2]), lim1.logg[ul_g1], lim2.logg[ul_g2], r1, r2, r1, r2)
    end
end
function (table::YBCTable)(Teff::Real, logg::Real, Mdot::Real)
    limits = table.limits
    tables = table.tables
    arg = (Teff=Teff, logg=logg, Mdot=Mdot) # argument NamedTuple to pass into the individual tables
    flags = sacollect(SVector{length(limits), Bool}, check_vals(lim, Teff, logg, Mdot) for lim in limits)
    c = count(flags)
    if c == 1
        idx = findfirst(flags)
        println(idx)
        return _get_table(tables, idx)(arg)
        # t, _ =  _get_table(tables, limits, idx)
        # return t(arg)
    elseif c == 2
        idxs = findall(flags)
        println(idxs)
        i, j = idxs[1], idxs[2]
        t1 = _get_table(tables, i)
        t2 = _get_table(tables, j)
        # Call interpolation between the two libraries
        return _itp_tables(t1, t2, limits[i], limits[j], Teff, logg, Mdot)
    elseif c > 2
        @warn "Provided combination of Teff=$Teff, logg=$logg, Mdot=$Mdot matched more than 2 YBC sub-libraries. Smoothing interpolation between libraries is only supported for up to 2."
    else # c == 0
        error("Provided Teff=$Teff, logg=$logg, Mdot=$Mdot did not match a YBC sub-library.")
    end
    # return flags
end
(table::YBCTable)(arg) = table(_parse_teff(arg), _parse_logg(arg), _parse_Mdot(arg))
# (table::WMbasicYBCTable)(model::Bjorklund2021MassLoss, arg) = table(_parse_teff(arg), _parse_logg(arg), _parse_Mdot(arg, Z(table), model))
# # Data are naturally Float32 -- convert hardware numeric args for faster evaluation and guarantee Float32 output
(table::YBCTable)(Teff::HardwareNumeric, logg::HardwareNumeric, Mdot::HardwareNumeric) = table(convert(dtype, Teff), convert(dtype, logg), convert(dtype, Mdot))
# # Methods to fix method ambiguities
# (::YBCTable)(::AbstractArray{<:Real}) = throw(ArgumentError("Requires at least 2 input arrays (Teff, logg)."))
# (::YBCTable)(::Type{Table}) = throw(ArgumentError("Requires at least 2 input arrays (Teff, logg)."))
# to broadcast over both teff and logg, you do table.(teff, logg')


 export YBCGrid, YBCTable

##########################################################################

end # Module