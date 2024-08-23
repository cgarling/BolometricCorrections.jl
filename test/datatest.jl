import FITSIO: FITS, read_header
import TypedTables: Table, columnnames

function read_fitstable(filename)
    result = FITS(filename) do f
        Table(f[2])
    end
    return result
end

const dirpath = "/home/cgarling/Work/Resources/stellar_models/YBC_tables/YBC/acs_wfc_202101/regrid"
dpath(filename) = joinpath(dirpath, filename)



# PHOENIX, naming convention ...Settl_M-2.0_a+0.4 has [M/H] = -2 and [alpha/H] = +0.4
a = read_fitstable(dpath("Avodonnell94Rv3.1BT-Settl_M-2.0_a+0.4.BC.fits"))
# CK2004 atlas9 model; readme.md says that `fm05k2odfnew` corresponds to
# `m05` means [M/H]=-0.05 (`p05` means [M/H]=0.05), but Avodonnell94Rv3.1YBC.list
# indicates that `m05` means [M/H] = -0.5
b = read_fitstable(dpath("Avodonnell94Rv3.1fm25k2odfnew.BC.fits"))
# ATLAS12 models for 47 Tuc (two populations)
c = read_fitstable(dpath("Avodonnell94Rv3.147TucP1.BC.fits"))
# COMARCS C-stars; http://stev.oapd.inaf.it/atm/lrspe.html seems to give description
# of filename:
# o  = log(eps(O)/eps(H)) + 12*
# c  = log(eps(C)/eps(H)) + 12*
# fe = log(eps(Fe)/eps(H)) + 12*
# Lines in Avodonnell94Rv3.1YBC.list seem to be filename, Z, [Fe/H], C/O, <something> 
d1 = read_fitstable(dpath("Avodonnell94Rv3.1COMARCS_Cstars_o83029_c83241_fe70630.BC.fits"))
# COMARCS M-stars, I think `feh-2.00` actually means [M/H] = -2,
# as I believe the metallicity of the grids is defined in Z
# See StarFormationHistories.Z_from_MH)
d2 = read_fitstable(dpath("Avodonnell94Rv3.1COMARCS_Mstars_mstar_feh-2.00_bt2.dat.BC.fits"))
# PoWR models for WR stars
e = read_fitstable(dpath("Avodonnell94Rv3.1subsmc-wnegrid.200-80000.extended.ori.BC.fits"))
# WMbasic models
f = read_fitstable(dpath("Avodonnell94Rv3.1WM_Z0.0001Mdot5.BC.fits"))
# Tlusty models, not mentioned in Chen2019
g = read_fitstable(dpath("Avodonnell94Rv3.1Tlusty_H0.2.BC.fits"))
