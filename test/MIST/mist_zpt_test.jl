using BolometricCorrections
using TypedTables: Table
using Test

const zpt = BolometricCorrections.MIST.zpt # Constant instance of BolometricCorrections.MIST.MISTZeropoints


@test Table(zpt) == zpt.table
@test zpt("WISE_W1") isa NamedTuple
@test zpt(@view "TWISE_W1"[2:end]) isa NamedTuple # Test with SubString
@test zpt.WISE_W1 == zpt("WISE_W1")
@test filternames(zpt) isa AbstractVector{<:AbstractString}
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

    # Bolometric corrections are linear in Av between the nodes returned by Av_nodes
    @test Av_nodes(grid) == Av_nodes(MISTv1BCGrid) == collect(BolometricCorrections.MIST.gridinfov1.Av)
    t4, t6 = grid(-1.01, 0.4), grid(-1.01, 0.6)
    @test grid(-1.01, 0.55)(5000.0, 4.0) ≈ t4(5000.0, 4.0) .+ 0.75 .* (t6(5000.0, 4.0) .- t4(5000.0, 4.0))

    # Grids with a subset of the filters give the same BCs for those filters, in the order given
    sub = select_filters(grid, ("F200W", :F090W))
    @test sub isa MISTv1BCGrid
    @test filternames(sub) == filternames(sub(-1.01, 0.43)) == (:F200W, :F090W)
    i = [findfirst(==(f), filternames(grid)) for f in (:F200W, :F090W)]
    @test sub(-1.01, 0.43)([3000.0, 5000.0], [1.0, 4.0]) == grid(-1.01, 0.43)([3000.0, 5000.0], [1.0, 4.0])[i, :]
    @test_throws "not found" select_filters(grid, ("F200W", "F999W"))
    @test_throws ArgumentError select_filters(grid, ("F200W", "F200W"))

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

    @test Av_nodes(grid) == Av_nodes(MISTv2BCGrid) == collect(BolometricCorrections.MIST.gridinfov2.Av)
    t4, t6 = grid(-1.01, 0.0, 0.4), grid(-1.01, 0.0, 0.6)
    @test grid(-1.01, 0.0, 0.55)(5000.0, 4.0) ≈ t4(5000.0, 4.0) .+ 0.75 .* (t6(5000.0, 4.0) .- t4(5000.0, 4.0))

    fs = filternames(grid)[[3, 1]]
    sub = select_filters(grid, fs)
    @test sub isa MISTv2BCGrid
    @test filternames(sub) == filternames(sub(-1.01, 0.2, 0.43)) == fs
    @test sub(-1.01, 0.2, 0.43)([3000.0, 5000.0], [1.0, 4.0]) == grid(-1.01, 0.2, 0.43)([3000.0, 5000.0], [1.0, 4.0])[[3, 1], :]

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