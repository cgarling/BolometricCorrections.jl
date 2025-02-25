using BolometricCorrections.MIST: MISTChemistry, X, X_phot, Y, Y_phot, Z, Z_phot, MH
using Test: @test

# Presently the MIST chemistry is a little wonky because the Asplund protostellar abundances
# (last row in Table 4 of Asplund 2009) have X + Y + Z = 0.9999, not 1
# Fixing this by dividing the table values by 0.9999 to renormalize sum X + Y + Z = 1.

# const rtol = 1e-3
const chem = MISTChemistry()
@test X(chem) + Y(chem) + Z(chem) ≈ 1
@test X_phot(chem) + Y_phot(chem) + Z_phot(chem) ≈ 1

@test X(chem, Z(chem)) ≈ X(chem) # rtol=rtol
@test X_phot(chem, Z_phot(chem)) ≈ X_phot(chem) # rtol=rtol

@test Y(chem, Z(chem)) ≈ Y(chem)
@test Y_phot(chem, Z_phot(chem)) ≈ Y_phot(chem)

@test MH(chem, 1e-4) ≈ -2.173234648003277 # rtol=rtol
@test MH(chem, Z(chem, 0.5)) ≈ 0.5 # rtol=rtol
@test Z(chem, MH(chem, 1e-2)) ≈ 1e-2 # rtol=rtol
