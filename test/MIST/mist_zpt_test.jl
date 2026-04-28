using BolometricCorrections
using TypedTables: Table
using Test

const zpt = BolometricCorrections.MIST.zpt # Constant instance of BolometricCorrections.MIST.MISTZeropoints


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

@testset "v1" begin
    grid = MISTv1BCGrid("JWST")
    table1 = grid(-1.01, 0.01)
    @test zeropoints(grid) === zpt
    @test zeropoints(table1) === zpt

    @test gridname(MISTv1BCGrid) isa String
    @test gridname(grid) isa String
    @test gridname(MISTv1BCTable) isa String
    @test gridname(table1) isa String

    for feh in range(extrema(BolometricCorrections.MIST.gridinfov1.feh)...; step=0.1)
        for Av in range(extrema(BolometricCorrections.MIST.gridinfov1.Av)...; step=0.1)
            table = grid(feh, Av)
            @test table isa MISTv1BCTable
            @test table(3000, 1.0) isa AbstractVector # technically SVector but ...
        end
    end
end

@testset "v2" begin
    grid = MISTv2BCGrid("JWST")
    table1 = grid(-1.01, 0.0, 0.01)
    @test zeropoints(grid) === zpt
    @test zeropoints(table1) === zpt

    @test gridname(MISTv2BCGrid) isa String
    @test gridname(grid) isa String
    @test gridname(MISTv2BCTable) isa String
    @test gridname(table1) isa String

    for feh in range(extrema(BolometricCorrections.MIST.gridinfov2.feh)...; step=0.1)
        for Av in range(extrema(BolometricCorrections.MIST.gridinfov2.Av)...; step=0.1)
            for afe in range(extrema(BolometricCorrections.MIST.gridinfov2.afe)...; step=0.1)
                table = grid(feh, afe, Av)
                @test table isa MISTv2BCTable
                @test table(3000, 1.0) isa AbstractVector # technically SVector but ...
            end
        end
    end
end