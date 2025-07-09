module YBC

# using ..BolometricCorrections: @compat
using Compat: @compat
import CSV
using TypedTables: Table

# Code to initialize data storage mechanisms
include("init.jl")

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

include("PHOENIX.jl")
using .PHOENIX
@compat public PHOENIX


end # Module