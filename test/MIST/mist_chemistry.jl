using BolometricCorrections.MIST: MISTChemistry, MISTBCGrid, X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry
using Test: @test

# Presently the MIST chemistry is a little wonky because the Asplund protostellar abundances
# (last row in Table 4 of Asplund 2009) have X + Y + Z = 0.9999, not 1
# Fixing this by dividing the table values by 0.9999 to renormalize sum X + Y + Z = 1.

# const rtol = 1e-3
const chem = MISTChemistry()
@test X(chem) + Y(chem) + Z(chem) ≈ 1
@test X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1

@test X(chem, Z(chem)) ≈ X(chem) # rtol=rtol

@test Y(chem, Z(chem)) ≈ Y(chem)

@test MH(chem, 1e-4) ≈ -2.173234648003277 # rtol=rtol
@test MH(chem, Z(chem, 0.5)) ≈ 0.5 # rtol=rtol
@test Z(chem, MH(chem, 1e-2)) ≈ 1e-2 # rtol=rtol

grid = MISTBCGrid("JWST")
for feh in range(extrema(grid).feh...; step=10)
    table = grid(feh, 0.0)
    @test MH(table) == feh
    tchem = chemistry(table)
    @test Z(table) ≈ Z(tchem, MH(table))
    @test Y(table) ≈ Y(tchem, Z(table))
    @test X(table) ≈ X(tchem, Z(table))
end
