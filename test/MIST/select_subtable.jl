using Test: @test
using TypedTables: Table
import BolometricCorrections.MIST

# Test that select_subtable correctly selects rows with feh == feh and Av == Av
table = MIST.MISTBCGrid("JWST")
for feh in MIST.gridinfo.feh
    for Av in MIST.gridinfo.Av
        subtable = MIST.select_subtable(Table(table), feh, Av)
        @test length(subtable) == length(MIST.gridinfo.Teff) * length(MIST.gridinfo.logg)
        @test all(@. subtable.feh == feh && subtable.Av == Av)
    end
end
