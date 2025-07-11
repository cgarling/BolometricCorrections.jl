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

grid = MISTBCGrid("JWST")
@test zeropoints(grid) === zpt
table = grid(-1.01, 0.01)
@test zeropoints(table) === zpt
for feh in range(extrema(BolometricCorrections.MIST.gridinfo.feh)...; step=0.1)
    for Av in range(extrema(BolometricCorrections.MIST.gridinfo.Av)...; step=0.1)
        table = grid(feh, Av)
        @test table isa MISTBCTable
        @test table(3000, 1.0) isa AbstractVector # technically SVector but ...
    end
end