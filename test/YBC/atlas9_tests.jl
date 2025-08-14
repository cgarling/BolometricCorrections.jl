using Test: @test
using BolometricCorrections
using StaticArrays: SVector

# Mixed types for mh and Av used to be incorrectly converted, so we test

grid = ATLAS9YBCGrid("acs_wfc")
table1 = grid(-1.0, 0f0)
# One regression test
@test table1(3500.0, 0.1) ≈ BolometricCorrections.YBC.dtype[-3.8445332, -2.9072084, -2.7502513, -1.5907421, -1.9208325, -1.2539452, -0.87817866, -0.4435694, -0.5165189, 0.12026592, 0.28986067, 0.6828744]
@test gridname(grid) isa String
@test gridname(ATLAS9YBCGrid) isa String
@test gridname(table1) isa String
@test gridname(ATLAS9YBCTable) isa String
for mh in Float32.(range(extrema(BolometricCorrections.YBC.ATLAS9.gridinfo.MH)...; step=0.1))
    for Av in Float64.(range(extrema(BolometricCorrections.YBC.ATLAS9.gridinfo.Av)...; step=0.1))
        table = grid(mh, Av)
        @test table isa ATLAS9YBCTable
        @test MH(table) ≈ mh
        chem = chemistry(table)
        @test Z(table) ≈ Z(chem, MH(table))
        @test Y(table) ≈ Y(chem, Z(table))
        @test X(table) ≈ X(chem, Z(table))
        @test table(4250.0, 1.0) isa SVector{12, BolometricCorrections.YBC.dtype}
    end
end