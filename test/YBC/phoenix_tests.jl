using Test: @test
using BolometricCorrections
using StaticArrays: SVector

# Mixed types for mh and Av used to be incorrectly converted, so we test

grid = PHOENIXYBCGrid("acs_wfc")
table1 = grid(-1.0, 0f0)
# One regression test
@test table1(3500.0, 0.1) ≈ BolometricCorrections.YBC.dtype[-3.6615107, -2.9102192, -2.8467753, -1.7102648, -2.0258598, -1.506531, -1.1820091, -0.9494957, -0.9117589, 0.14530413, 0.33153394, 0.73742735]
@test gridname(grid) isa String
@test gridname(PHOENIXYBCGrid) isa String
@test gridname(table1) isa String
@test gridname(PHOENIXYBCTable) isa String

# Test FeH / alphaFe at grid [M/H] nodes
_mh_afe = Dict(zip(BolometricCorrections.YBC.PHOENIX._mh, BolometricCorrections.YBC.PHOENIX._afe))

for (mh, expected_afe) in pairs(_mh_afe)
    table = grid(mh, 0.0f0)
    f_alpha = alpha_mass_fraction(chemistry(table))
    @test FeH(table) ≈ mh - log10(f_alpha * exp10(expected_afe) + (1 - f_alpha))
    @test alphaFe(table) ≈ expected_afe
end

for mh in Float32.(range(extrema(BolometricCorrections.YBC.PHOENIX.gridinfo.MH)...; step=0.1))
    for Av in Float64.(range(extrema(BolometricCorrections.YBC.PHOENIX.gridinfo.Av)...; step=0.1))
        table = grid(mh, Av)
        @test table isa PHOENIXYBCTable
        @test MH(table) ≈ mh
        chem = chemistry(table)
        @test Z(table) ≈ Z(chem, MH(table))
        @test Y(table) ≈ Y(chem, Z(table))
        @test X(table) ≈ X(chem, Z(table))
        @test table(3000.0, 1.0) isa SVector{12, BolometricCorrections.YBC.dtype}
    end
end