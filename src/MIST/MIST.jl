module MIST

# using ..BolometricCorrections: Table, columnnames # relative path for parent module
using CodecXz: XzDecompressorStream # Decompress downloaded BCs
import CSV
using DataDeps: register, DataDep, @datadep_str
import HDF5
using Printf: @sprintf # Formatted conversion of floats to strings
import Tar
import Tables # for Tables.matrix conversion
using TypedTables: Table, columnnames

# Initialization and datadeps
include("init.jl")


end # module
