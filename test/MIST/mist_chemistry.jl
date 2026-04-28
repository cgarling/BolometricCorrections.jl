using BolometricCorrections.MIST: MISTv1Chemistry, MISTv2Chemistry, MISTv1BCGrid, MISTv2BCGrid, X, X_phot, Y, Y_phot, Z, Z_phot, MH, FeH, alphaFe, chemistry, alpha_mass_fraction
using Test: @test, @testset

# Presently the MIST chemistry is a little wonky because the Asplund protostellar abundances
# (last row in Table 4 of Asplund 2009) have X + Y + Z = 0.9999, not 1
# Fixing this by dividing the table values by 0.9999 to renormalize sum X + Y + Z = 1.

# const rtol = 1e-3
@testset "MISTv1Chemistry" begin

    chem = MISTv1Chemistry()
    @test X(chem) + Y(chem) + Z(chem) ≈ 1
    @test X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1
    @test X(chem, Z(chem)) ≈ X(chem)
    @test Y(chem, Z(chem)) ≈ Y(chem)

    @test MH(chem, 1e-4) ≈ -2.173234648003277 # rtol=rtol
    @test MH(chem, Z(chem, 0.5)) ≈ 0.5 # rtol=rtol
    @test Z(chem, MH(chem, 1e-2)) ≈ 1e-2 # rtol=rtol

    grid = MISTv1BCGrid("JWST")
    for feh in range(extrema(grid).feh...; step=10)
        table = grid(feh, 0.0)
        @test MH(table) == feh
        @test FeH(table) == feh
        @test alphaFe(table) == 0
        tchem = chemistry(table)
        @test Z(table) ≈ Z(tchem, MH(table))
        @test Y(table) ≈ Y(tchem, Z(table))
        @test X(table) ≈ X(tchem, Z(table))
    end
end

@testset "MISTv2Chemistry" begin
    chem = MISTv2Chemistry()
    @test X(chem) + Y(chem) + Z(chem) ≈ 1
    @test X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1
    @test X(chem, Z(chem)) ≈ X(chem) # rtol=rtol
    @test Y(chem, Z(chem)) ≈ Y(chem)

    @test MH(chem, 1e-4) ≈ -2.292643973966403 # rtol=rtol
    @test MH(chem, Z(chem, 0.5)) ≈ 0.5 # rtol=rtol
    @test Z(chem, MH(chem, 1e-2)) ≈ 1e-2 # rtol=rtol

    grid = MISTv2BCGrid("JWST")
    for feh in range(extrema(grid).feh...; step=10)
        for afe in (-0.2, 0.0) # [α/Fe]
            table = grid(feh, afe, 0.0)
            @test FeH(table) == feh
            @test alphaFe(table) == afe
            f_alpha = alpha_mass_fraction(chemistry(table))
            @test MH(table) ≈ feh + log10(f_alpha * exp10(alphaFe(table)) + (1 - f_alpha))
            tchem = chemistry(table)
            @test Z(table) ≈ Z(tchem, MH(table))
            @test Y(table) ≈ Y(tchem, Z(table))
            @test X(table) ≈ X(tchem, Z(table))
        end
    end
end
