module YBC

# import ..MIST: mist_function # relative path for a “sibling” module
using ..BolometricCorrections: Table, columnnames # relative path for parent module

""""
    Calculates the abundance ratio ``C/O = N_C /N_O``, where ``N_X`` is the number density of element ``X``, from values `o` and `c` from COMARCS model filenames. Three argument version also returns [Fe/H].

COMARCS defines `o` and `c` as ``o = log(eps(O)/eps(H)) + 12``, ``c = log(eps(C)/eps(H)) + 12``, ``fe = log(eps(Fe)/eps(H)) + 12``, see [this page](http://stev.oapd.inaf.it/atm/lrspe.html). 
"""
function comarcs_CO(o, c)
    # ``[C/O] = log (N_C/N_O) − log (N_C/N_O)_{sol}``
    n_o = exp10(o - 12) # n_o / n_h
    n_c = exp10(c - 12) # n_c / n_h
    c_o = n_c / n_o
end
comarcs_CO(o, c, fe) = comarcs_CO(o, c), fe - 12
function comarcs_CO(fname::AbstractString)
    # Decimal place is after first digit and string descriptors are all 5 characters, e.g.
    # ..._o83029_c83241_fe70630.BC.fits, so we can just divide by 1e4
    o = parse(Float64, split(split(fname, "_o")[2], '_')[1]) / 1e4
    c = parse(Float64, split(split(fname, "_c")[2], '_')[1]) / 1e4
    fe = parse(Float64, split(split(fname, "_fe")[2], '.')[1]) / 1e4
    return comarcs_CO(o, c, fe)
end

# const test2 = "asdfg" # This is exported at the submodule level
# export test2
end # module
