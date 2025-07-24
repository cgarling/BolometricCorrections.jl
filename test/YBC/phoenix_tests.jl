using Test: @test
using BolometricCorrections
using StaticArrays: SVector

grid = PHOENIXYBCGrid("acs_wfc")
for mh in range(extrema(BolometricCorrections.YBC.PHOENIX.gridinfo.MH)...; step=0.1)
    for Av in range(extrema(BolometricCorrections.YBC.PHOENIX.gridinfo.Av)...; step=0.1)
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