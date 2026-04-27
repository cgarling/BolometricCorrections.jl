using Test: @test
using TypedTables: Table
import BolometricCorrections.MIST

@testset "v1" begin
    # Test that select_subtable correctly selects rows with feh == feh and Av == Av
    table = MIST.MISTBCGridv1("JWST")
    for feh in MIST.gridinfov1.feh
        for Av in MIST.gridinfov1.Av
            subtable = MIST.select_subtable_v1(Table(table), feh, Av)
            @test length(subtable) == length(MIST.gridinfov1.Teff) * length(MIST.gridinfov1.logg)
            # @test all(@. subtable.feh == feh && subtable.Av == Av) # Not working on julia 1.6
            @test all(@. (subtable.feh == feh) & (subtable.Av == Av))
        end
    end
end

@testset "v2" begin
    table = MIST.MISTBCGridv2("JWST")
    for feh in MIST.gridinfov2.feh
        for Av in MIST.gridinfov2.Av
            for afe in MIST.gridinfov2.afe
                subtable = MIST.select_subtable_v2(Table(table), feh, afe, Av)
                @test length(subtable) == length(MIST.gridinfov2.lgTef) * length(MIST.gridinfov2.logg)
                # @test all(@. subtable.feh == feh && subtable.Av == Av) # Not working on julia 1.6
                @test all(@. (subtable.feh == feh) & (subtable.Av == Av) & (subtable.afe == afe))
            end
        end
    end
end
