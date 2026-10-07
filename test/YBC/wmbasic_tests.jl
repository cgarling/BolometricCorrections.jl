using Test: @test
using BolometricCorrections
using BolometricCorrections.YBC.WMbasic: WMbasicYBCGrid, WMbasicYBCTable, gridinfo
using StaticArrays: SVector

# Mixed types for mh and Av used to be incorrectly converted, so we test

grid = WMbasicYBCGrid("acs_wfc")
table1 = grid(-1.0, 0f0)
@test table1(24_000.0, 3.14, 1e-6) ≈ BolometricCorrections.YBC.dtype[-2.1013548, -2.2280855, -2.3163576, -2.3767126, -2.333343, -2.3797524, -2.4334195, -2.3151662, -2.4037402, -2.5517056, -2.5801084, -2.671912]
@test gridname(WMbasicYBCGrid) isa String
@test gridname(grid) isa String
@test gridname(WMbasicYBCTable) isa String
@test gridname(table1) isa String
for mh in Float32.(range(extrema(grid).MH[1], extrema(grid).MH[2]; length=10))
    for Av in Float64.(range(extrema(grid).Av...; length=10))
        table = grid(mh, Av)
        @test table isa WMbasicYBCTable
        @test MH(table) ≈ mh
        @test FeH(table) == MH(table)
        @test alphaFe(table) == 0
        chem = chemistry(table)
        @test Z(table) ≈ Z(chem, MH(table))
        @test Y(table) ≈ Y(chem, Z(table))
        @test X(table) ≈ X(chem, Z(table))
        for Mdot in range(extrema(table).Mdot...; length=10)
            @test table(20_250.0, 3.51, Mdot) isa SVector{12, BolometricCorrections.YBC.dtype}
        end
    end
end
# At grid nodes, the interpolated BCs must reproduce the YBC FITS files in every filter. The filter and Mdot dimensions
# of the interpolation grid were once transposed, scrambling the BCs without changing any array shapes.
using FITSIO: FITS
using BolometricCorrections.YBC: pull_table
using BolometricCorrections.YBC.WMbasic: _wmbasic_filename
path = pull_table("acs_wfc", "YBC")
for (iZ, iMdot, Av) in ((6, 1, 0), (3, 2, 1), (4, 3, 5))
    Zval, Mdot = gridinfo.Z[iZ], gridinfo.Mdot[iMdot]
    table = grid(gridinfo.MH[iZ], Av)
    FITS(joinpath(path, "regrid", _wmbasic_filename(Zval, Mdot))) do f
        logTeff, logg = read(f[2], "logTeff"), read(f[2], "logg")
        suffix = iszero(Av) ? "" : "_Av" * string(Av)
        for i in (1, length(logg) ÷ 2, length(logg))
            @test table(exp10(logTeff[i]), logg[i], Mdot) ≈ [read(f[2], String(filt) * suffix)[i] for filt in filternames(grid)] rtol=1e-4
        end
    end
end
