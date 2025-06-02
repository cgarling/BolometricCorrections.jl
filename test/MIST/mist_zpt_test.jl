using BolometricCorrections
zpt = BolometricCorrections.MIST.zpt # Constant instance of BolometricCorrections.MIST.MISTZeropoints
using TypedTables: Table

using Test

@test Table(zpt) == zpt.table
@test zpt("WISE_W1") isa NamedTuple
@test zpt(@view "TWISE_W1"[2:end]) isa NamedTuple # Test with SubString
@test zpt.WISE_W1 == zpt("WISE_W1")
@test filternames(zpt) isa Vector{<:AbstractString}
# WISE_W1 is in Vega system; test conversions
@test vegamags(zpt, "WISE_W1", 1.0) == 1.0
@test abmags(zpt, "WISE_W1", 1.0) == 1.0 + zpt("WISE_W1").VegaAB
@test stmags(zpt, "WISE_W1", 1.0) == 1.0 + zpt("WISE_W1").VegaST
@test Mbol(zpt) isa Number
@test Lbol(zpt) isa Number
