using Test: @test
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid, KoesterWDYBCTable, gridinfo, chemistry, dtype
using StaticArrays: SVector

@test ismissing(chemistry(KoesterWDYBCGrid))
@test ismissing(chemistry(KoesterWDYBCTable))

grid = KoesterWDYBCGrid("acs_wfc")
for Av in range(extrema(gridinfo.Av)...; step=0.1)
    table = grid(Av)
    @test table isa KoesterWDYBCTable
    @test table(10_250.0, 7.33) isa SVector{12, dtype}
end
