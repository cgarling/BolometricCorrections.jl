using Test: @test
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid, KoesterWDYBCTable, gridinfo, chemistry, gridname, dtype
using StaticArrays: SVector

@test ismissing(chemistry(KoesterWDYBCGrid))
@test ismissing(chemistry(KoesterWDYBCTable))

grid = KoesterWDYBCGrid("acs_wfc")
table1 = grid(0.0)
@test gridname(grid) isa String
@test gridname(KoesterWDYBCGrid) isa String
@test gridname(table1) isa String
@test gridname(KoesterWDYBCTable) isa String
for Av in range(extrema(gridinfo.Av)...; step=0.1)
    table = grid(Av)
    @test table isa KoesterWDYBCTable
    @test table(10_250.0, 7.33) isa SVector{12, dtype}
end
