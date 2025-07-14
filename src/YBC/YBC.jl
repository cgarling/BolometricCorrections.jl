module YBC

# using ..BolometricCorrections: @compat
using ..BolometricCorrections: AbstractChemicalMixture
import ..BolometricCorrections: X, X_phot, Y, Y_phot, Y_p, Z, Z_phot, MH, chemistry

using Compat: @compat
import CSV
using TypedTables: Table

# Code to initialize data storage mechanisms
include("init.jl")

"""YBC data are stored in FITS files with natural datatype Float32."""
const dtype = Float32

"""
    parse_filterinfo(f::AbstractString)
Given a `filter.info` file from YBC, return a `TypedTables.Table` with the corresponding information. The columns are described as follows by YBC:

`filter.info`: infomation for the filters. The columns are filter number, filter name, filter file used, effective wavelength (in Angstrom), band width, absolute flux of Vega (in unit of erg cm-2 s-1 A-1), photometric system (AB, Vega, ST or TG), type of device (energy counting or number counting for CCD), reference magnitude of Vega applied, irrelevant value, comments.
"""
function parse_filterinfo(f::AbstractString)
    # Contains final column with # <comments>, we don't want these so filter first
    dtype = Float64
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
This function expects the available [M/H] values to be `gridinfo.MH` and A_v values to be `gridinfo.Av`.
Submodules define 2-argument `check_vals` that use their own `gridinfo` transparently.
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
julia> using StellarTracks.PARSEC: PARSECChemistry, X, Y, Z, X_phot, Y_phot, Z_phot, MH;

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

end # Module