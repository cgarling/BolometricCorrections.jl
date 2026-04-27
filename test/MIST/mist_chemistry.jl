using BolometricCorrections.MIST: MISTChemistryv1, MISTChemistryv2, MISTBCGridv1, MISTBCGridv2, X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry
using Test: @test, @testset

# Presently the MIST chemistry is a little wonky because the Asplund protostellar abundances
# (last row in Table 4 of Asplund 2009) have X + Y + Z = 0.9999, not 1
# Fixing this by dividing the table values by 0.9999 to renormalize sum X + Y + Z = 1.

# const rtol = 1e-3
@testset "MISTChemistryv1" begin

    chem = MISTChemistryv1()
    @test X(chem) + Y(chem) + Z(chem) ≈ 1
    @test X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1
    @test X(chem, Z(chem)) ≈ X(chem)
    @test Y(chem, Z(chem)) ≈ Y(chem)

    @test MH(chem, 1e-4) ≈ -2.173234648003277 # rtol=rtol
    @test MH(chem, Z(chem, 0.5)) ≈ 0.5 # rtol=rtol
    @test Z(chem, MH(chem, 1e-2)) ≈ 1e-2 # rtol=rtol

    grid = MISTBCGridv1("JWST")
    for feh in range(extrema(grid).feh...; step=10)
        table = grid(feh, 0.0)
        @test MH(table) == feh
        tchem = chemistry(table)
        @test Z(table) ≈ Z(tchem, MH(table))
        @test Y(table) ≈ Y(tchem, Z(table))
        @test X(table) ≈ X(tchem, Z(table))
    end
end

@testset "MISTChemistryv2" begin
    chem = MISTChemistryv2()
    @test X(chem) + Y(chem) + Z(chem) ≈ 1
    @test X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1
    @test X(chem, Z(chem)) ≈ X(chem) # rtol=rtol
    @test Y(chem, Z(chem)) ≈ Y(chem)

    @test MH(chem, 1e-4) ≈ -2.292643973966403 # rtol=rtol
    @test MH(chem, Z(chem, 0.5)) ≈ 0.5 # rtol=rtol
    @test Z(chem, MH(chem, 1e-2)) ≈ 1e-2 # rtol=rtol

    grid = MISTBCGridv2("JWST")
    for feh in range(extrema(grid).feh...; step=10)
        for afe in (-0.2, 0.0) # [α/Fe]
            table = grid(feh, afe, 0.0)
            @test MH(table) == feh
            tchem = chemistry(table)
            @test Z(table) ≈ Z(tchem, MH(table))
            @test Y(table) ≈ Y(tchem, Z(table))
            @test X(table) ≈ X(tchem, Z(table))
        end
    end
end
