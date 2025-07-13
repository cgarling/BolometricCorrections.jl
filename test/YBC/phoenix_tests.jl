using Test: @test
using BolometricCorrections

grid = PHOENIXYBCGrid("acs_wfc")
for mh in range(extrema(BolometricCorrections.YBC.PHOENIX.gridinfo.MH)...; step=0.1)
    for Av in range(extrema(BolometricCorrections.YBC.PHOENIX.gridinfo.Av)...; step=0.1)
        table = grid(mh, Av)
        @test table isa PHOENIXYBCTable
        @test table(3000.0, 1.0) isa AbstractVector # technically SVector but ...
    end
end