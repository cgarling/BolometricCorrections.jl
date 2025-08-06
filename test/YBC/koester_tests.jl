using Test: @test
using BolometricCorrections.YBC.KoesterWD: KoesterWDYBCGrid, KoesterWDYBCTable, gridinfo, chemistry, gridname, dtype
using StaticArrays: SVector

@test ismissing(chemistry(KoesterWDYBCGrid))
@test ismissing(chemistry(KoesterWDYBCTable))

grid = KoesterWDYBCGrid("acs_wfc")
table1 = grid(0.0)
# One regression test
@test table1(6000.0, 7.2) ≈ dtype[-0.6667951, -0.46817672, -0.3289503, -0.07198871, -0.16654676, 0.024434362, 0.17963289, 0.33361864, 0.32127845, 0.522821, 0.5634577, 0.6676144]
@test gridname(grid) isa String
@test gridname(KoesterWDYBCGrid) isa String
@test gridname(table1) isa String
@test gridname(KoesterWDYBCTable) isa String
for Av in range(extrema(gridinfo.Av)...; step=0.1)
    table = grid(Av)
    @test table isa KoesterWDYBCTable
    @test table(10_250.0, 7.33) isa SVector{12, dtype}
end
