using Test: @test
using BolometricCorrections
using BolometricCorrections.YBC.WMbasic: WMbasicYBCGrid, WMbasicYBCTable, gridinfo
using StaticArrays: SVector

# Mixed types for mh and Av used to be incorrectly converted, so we test

grid = WMbasicYBCGrid("acs_wfc")
table1 = grid(-1.0, 0f0)
@test table1(24_000.0, 3.14, 1e-6) ≈ BolometricCorrections.YBC.dtype[-2.208993, -2.3101966, -2.28842, -2.5616233, -2.2280855, -2.3333433, -2.3151662, -2.5801084, -2.2240276, -2.3268185, -2.2955875, -2.5536644]
@test gridname(WMbasicYBCGrid) isa String
@test gridname(grid) isa String
@test gridname(WMbasicYBCTable) isa String
@test gridname(table1) isa String
for mh in Float32.(range(extrema(grid).MH[1], extrema(grid).MH[2]; length=10))
    for Av in Float64.(range(extrema(grid).Av...; length=10))
        table = grid(mh, Av)
        @test table isa WMbasicYBCTable
        @test MH(table) ≈ mh
        chem = chemistry(table)
        @test Z(table) ≈ Z(chem, MH(table))
        @test Y(table) ≈ Y(chem, Z(table))
        @test X(table) ≈ X(chem, Z(table))
        for Mdot in range(extrema(table).Mdot...; length=10)
            @test table(20_250.0, 3.51, Mdot) isa SVector{12, BolometricCorrections.YBC.dtype}
        end
    end
end