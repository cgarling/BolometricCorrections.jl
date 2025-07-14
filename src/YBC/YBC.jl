module YBC

# using ..BolometricCorrections: @compat
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