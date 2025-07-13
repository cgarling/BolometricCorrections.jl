using Test: @test
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid, KoesterWDYBCTable, gridinfo, chemistry

@test ismissing(chemistry(KoesterWDYBCGrid))
@test ismissing(chemistry(KoesterWDYBCTable))

grid = KoesterWDYBCGrid("acs_wfc")
for Av in range(extrema(gridinfo.Av)...; step=0.1)
    table = grid(Av)
    @test table isa KoesterWDYBCTable
    @test table(10_250, 7.33) isa AbstractVector # technically SVector but ...
end
