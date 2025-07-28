using Test: @test
using BolometricCorrections
using StaticArrays: SVector

grid = ATLAS9YBCGrid("acs_wfc")
table1 = grid(-1.0, 0.0)
@test gridname(grid) isa String
@test gridname(ATLAS9YBCGrid) isa String
@test gridname(table1) isa String
@test gridname(ATLAS9YBCTable) isa String
for mh in range(extrema(BolometricCorrections.YBC.ATLAS9.gridinfo.MH)...; step=0.1)
    for Av in range(extrema(BolometricCorrections.YBC.ATLAS9.gridinfo.Av)...; step=0.1)
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