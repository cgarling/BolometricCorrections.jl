using Test: @test
using BolometricCorrections
using BolometricCorrections.YBC.WMbasic: WMbasicYBCGrid, WMbasicYBCTable, gridinfo

grid = WMbasicYBCGrid("acs_wfc")
for mh in range(extrema(grid).MH[1], extrema(grid).MH[2]; length=10)
    for Av in range(extrema(grid).Av...; length=10)
        table = grid(mh, Av)
        @test table isa WMbasicYBCTable
        @test MH(table) ≈ mh
        chem = chemistry(table)
        @test Z(table) ≈ Z(chem, MH(table))
        @test Y(table) ≈ Y(chem, Z(table))
        @test X(table) ≈ X(chem, Z(table))
        for Mdot in range(extrema(table).Mdot...; length=10)
            @test table(20_250.0, 3.51, Mdot) isa AbstractVector # technically SVector but ...
        end
    end
end