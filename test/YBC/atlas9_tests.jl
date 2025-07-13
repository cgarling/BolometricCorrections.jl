using Test: @test
using BolometricCorrections

grid = ATLAS9YBCGrid("acs_wfc")
for mh in range(extrema(BolometricCorrections.YBC.ATLAS9.gridinfo.MH)...; step=0.1)
    for Av in range(extrema(BolometricCorrections.YBC.ATLAS9.gridinfo.Av)...; step=0.1)
        table = grid(mh, Av)
        @test table isa ATLAS9YBCTable
        @test table(4250.0, 1.0) isa AbstractVector # technically SVector but ...
    end
end