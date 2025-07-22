using Test: @test
using BolometricCorrections
using BolometricCorrections.YBC

"""
    is_between(a, b, c) = (min(b, c) ≤ a ≤ max(b, c)) || (a ≈ b ≈ c)
Determine whether `a` is between `b` and `c` without knowing whether `b` or `c` is larger.
When working with Float32, if `b` and `c` are approximately equal, then the interpolation
that creates `a` also be approximately equal, but may suffer from rounding error. To account for this
possibility, the second clause is added so this function returns `true` if all arguments
are approximately equal.
"""
is_between(a, b, c) = (min(b, c) ≤ a ≤ max(b, c)) || (a ≈ b ≈ c)

grid = YBCGrid("acs_wfc")
# for mh in range(extrema(grid).MH...; step=0.1)
for mh in range(-2, 0; step=0.1)
    for Av in range(extrema(grid).Av...; step=0.1)
        # println(mh, " ", Av)
        table = grid(mh, Av)
        @test table isa YBCTable

        #########################
        #########################
        # Test WMbasic -- atlas9
        wmbasic_transitions = table.transitions.wmbasic
        #########################
        # Tests for values on boundaries
        # For Teff, Mdot on lower bounds, should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1], 4.15, wmbasic_transitions.Mdot[1]
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For Teff, Mdot on upper bounds, should return wmbasic result
        Teff, logg, Mdot = wmbasic_transitions.Teff[2], 4.15, wmbasic_transitions.Mdot[2]
        @test table(Teff, logg, Mdot) ≈ table.tables.wmbasic(Teff, logg, Mdot)
        # For Teff on lower bound, Mdot on upper bound, should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1], 4.15, wmbasic_transitions.Mdot[2]
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For Teff on upper bound, Mdot on lower bound, should return interpolation between atlas9, wmbasic
        Teff, logg, Mdot = wmbasic_transitions.Teff[2], 4.15, wmbasic_transitions.Mdot[1]
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.atlas9(Teff, logg), table.tables.wmbasic(Teff, logg, Mdot)
        @test mapreduce(is_between, &, r1, r2, r3) 

        #########################
        # Tests within interp region
        # For lower_bound < Teff < upper_bound, Mdot on lower bound, should interpolate between atlas9, wmbasic
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] + 100, 4.15, wmbasic_transitions.Mdot[1]
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.atlas9(Teff, logg), table.tables.wmbasic(Teff, logg, Mdot)
        @test mapreduce(is_between, &, r1, r2, r3) # Test element-wise that each entry of table(Teff, logg, Mdot) is between the atlas9 and wmbasic result
        # For lower_bound < Teff < upper_bound, Mdot on upper bound, should interpolate between atlas9, wmbasic
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] + 100, 4.15, wmbasic_transitions.Mdot[2]
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.atlas9(Teff, logg), table.tables.wmbasic(Teff, logg, Mdot)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For Teff on lower bound, lower bound < Mdot < upper bound, should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1], 4.15, wmbasic_transitions.Mdot[1] * 1.2
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For Teff on upper bound, lower bound < Mdot < upper bound, should interpolate between atlas9, wmbasic
        Teff, logg, Mdot = wmbasic_transitions.Teff[2], 4.15, wmbasic_transitions.Mdot[1] * 1.2
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.atlas9(Teff, logg), table.tables.wmbasic(Teff, logg, Mdot)
        @test mapreduce(is_between, &, r1, r2, r3) # Test element-wise that each entry of table(Teff, logg, Mdot) is between the atlas9 and wmbasic result
        # For lower_bound < Teff < upper_bound *and* lower bound < Mdot < upper bound, should interpolate between atlas9, wmbasic
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] + 100, 4.15, wmbasic_transitions.Mdot[1] * 1.2
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.atlas9(Teff, logg), table.tables.wmbasic(Teff, logg, Mdot)
        @test mapreduce(is_between, &, r1, r2, r3) # Test element-wise that each entry of table(Teff, logg, Mdot) is between the atlas9 and wmbasic result

        #########################
        # Tests beyond interp region
        # For Mdot > upper bound, lower bound < Teff < upper bound, should interpolate between atlas9, wmbasic
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] + 100, 4.15, wmbasic_transitions.Mdot[2] * 1.2
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.atlas9(Teff, logg), table.tables.wmbasic(Teff, logg, Mdot)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For Mdot < lower bound, lower bound < Teff < upper bound, should return atlas9
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] + 100, 4.15, wmbasic_transitions.Mdot[1] * 0.8
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For Mdot > upper bound, Teff < lower bound, should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] - 100, 4.15, wmbasic_transitions.Mdot[2] * 1.2
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For Mdot = upper bound, Teff < lower bound,  should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] - 100, 4.15, wmbasic_transitions.Mdot[2]
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # These two tests technically pass into the phoenix - atlas regime because Mdot < lower bound,
        # but wmbasic_transitions.Teff[1] is outside the phoenix regime so it always returns atlas9
        # For Mdot < lower bound, Teff > lower bound should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] + 100, 4.15, wmbasic_transitions.Mdot[1] * 0.8
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For Mdot < lower bound, Teff < lower bound should return atlas9 result
        Teff, logg, Mdot = wmbasic_transitions.Teff[1] - 100, 4.15, wmbasic_transitions.Mdot[1] * 0.8
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)


        #########################
        #########################
        # Test Koester -- PHOENIX
        koester_transitions = table.transitions.koester
        #########################
        # Tests for values on boundaries
        # For Teff, logg on lower bounds, should return phoenix result
        Teff, logg, Mdot = koester_transitions.Teff[1], koester_transitions.logg[1], 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # For Teff, logg on upper bounds, should return koester result
        Teff, logg, Mdot = koester_transitions.Teff[2], koester_transitions.logg[2], 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.koester(Teff, logg)
        # For Teff on lower bound, logg on upper bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[1], koester_transitions.logg[2], 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For Teff on upper bound, logg on lower bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[2], koester_transitions.logg[1], 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)

        #########################
        # Tests within interp region
        # For lower_bound < Teff < upper_bound, logg on lower bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[1] + 100, koester_transitions.logg[1], 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For lower_bound < Teff < upper_bound, logg on upper bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[1] + 100, koester_transitions.logg[2], 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For Teff on lower bound, lower bound < logg < upper bound, should return phoenix result
        Teff, logg, Mdot = koester_transitions.Teff[1], koester_transitions.logg[1] * 1.2, 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # For Teff on upper bound, lower bound < logg < upper bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[2], koester_transitions.logg[1] * 1.2, 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For lower_bound < Teff < upper_bound *and* lower bound < logg < upper bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[1] + 100, koester_transitions.logg[1] * 1.2, 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)

        #########################
        # Tests beyond interp region
        # For logg > upper bound, lower bound < Teff < upper bound, should interpolate between phoenix, koester
        Teff, logg, Mdot = koester_transitions.Teff[1] + 100, koester_transitions.logg[2] * 1.2, 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.koester(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)
        # For logg > upper bound, Teff < lower bound, should return phoenix result
        Teff, logg, Mdot = koester_transitions.Teff[1] - 100, koester_transitions.logg[2] * 1.2, 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # For logg = upper bound, Teff < lower bound,  should return phoenix result
        Teff, logg, Mdot = koester_transitions.Teff[1] - 100, koester_transitions.logg[2], 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # When logg < koester_transitions.logg[1] and Teff < koester_transitions.Teff[2],
        # we are in the atlas -- phoenix regime; test this in next block ...

        # For Mdot < lower bound, Teff > lower bound should return phoenix result
        # Teff, logg, Mdot = koester_transitions.Teff[1] + 100, 4.15, koester_transitions.Mdot[1] * 0.8
        # @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # # For Mdot < lower bound, Teff < lower bound should return phoenix result
        # Teff, logg, Mdot = koester_transitions.Teff[1] - 100, 4.15, koester_transitions.Mdot[1] * 0.8
        # @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # For logg < lower bound, lower bound < Teff < upper bound, should return phoenix
        # Note this only works because koester_transitions.Teff[1] + 100 does not reach the phoenix -- atlas transition region
        # Teff, logg, Mdot = koester_transitions.Teff[1] + 100, koester_transitions.logg[1] * 0.8, 0.0
        # @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)


        #########################
        #########################
        # Test PHOENIX -- ATLAS9
        phoenix_transitions = table.transitions.phoenix
        #########################
        # Tests for values on boundaries
        # For Teff on lower bounds, should return phoenix result
        Teff, logg, Mdot = phoenix_transitions.Teff[1], 2.0, 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.phoenix(Teff, logg)
        # For Teff on upper bounds, should return atlas9 result
        Teff, logg, Mdot = phoenix_transitions.Teff[2], 2.0, 0.0
        @test table(Teff, logg, Mdot) ≈ table.tables.atlas9(Teff, logg)
        # For lower bound < Teff < upper bound, should interpolate beween PHOENIX, ATLAS9
        Teff, logg, Mdot = phoenix_transitions.Teff[1] + 100, 2.0, 0.0
        r1, r2, r3  = table(Teff, logg, Mdot), table.tables.phoenix(Teff, logg), table.tables.atlas9(Teff, logg)
        @test mapreduce(is_between, &, r1, r2, r3)
    end
end