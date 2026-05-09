# Derivation of the solar alpha-element mass fraction for PARSECChemistry.
#
# The YBC bolometric corrections for PARSEC/PADOVA isochrones use the chemical mixture
# defined in §4 of:
#   Bressan, A. et al. (2012), MNRAS, 427, 127–145
#   DOI: 10.1111/j.1365-2966.2012.21948.x
#   ADS: https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract
#
# Bressan+2012 adopt Grevesse & Sauval (1998) [GS98] as the base solar mixture but replace
# selected element abundances with values from Caffau et al. and co-author papers, as listed
# explicitly in Table 1 of Bressan+2012:
#   C=8.50  Caffau et al. (2010)
#   N=7.86  Caffau et al. (2009)
#   O=8.76  Caffau et al. (2008b)
#   P=5.46  Caffau et al. (2007b)
#   S=7.16  Caffau & Ludwig (2007)
#   Fe=7.52 Caffau et al. (2011)
# All other elements (including Ne, Mg, Si, Ar, Ca, Ti) come from GS98.
# Caffau et al. overview paper:
#   Caffau, E. et al. (2011), Solar Physics, 268, 255–269
#   DOI: 10.1007/s11207-010-9541-4
#   ADS: https://ui.adsabs.harvard.edu/abs/2011SoPh..268..255C/abstract
#
# The alpha-element mass fraction f_α is defined as the fraction of total metal (Z) mass
# contributed by alpha elements: O, Ne, Mg, Si, S, Ar, Ca, Ti.
# It enters the Salaris et al. (1993) conversion between [M/H] and [Fe/H]:
#   [M/H] = [Fe/H] + log10(f_α · 10^[α/Fe] + (1 - f_α))
#
# For each element i, the mass contribution is proportional to A_i · 10^(ε_i - 12),
# where A_i is the atomic mass number and ε_i = log10(N_i / N_H) + 12 is the standard
# logarithmic abundance on the scale where log(N_H) + 12 = 12.
#
# Representative metals are included below, sufficient to evaluate f_α.
# Omitted trace elements contribute negligibly to Z and do not affect f_α at the
# precision retained here.

using Test: @test
using BolometricCorrections.YBC: PARSECChemistry, alpha_mass_fraction

# ────────────────────────────────────────────────────────────────────────────
# Bressan+2012 solar mixture, Table 1 of Bressan et al. (2012):
# GS98 base, with specific elements replaced by Caffau et al. sub-papers.
# Ne, Mg, Si, Ar, Ca, Ti are NOT in Table 1 and therefore come from GS98.
# Sources:
#   GS98       : Grevesse & Sauval (1998), Table 1
#   Caffau+2010 : C
#   Caffau+2009 : N
#   Caffau+2008b: O
#   Caffau+2007b: P
#   CaffauLudwig+2007: S
#   Caffau+2011 : Fe
# Columns: element => (atomic mass number A, log_epsilon)
# ────────────────────────────────────────────────────────────────────────────
const _bressan2012 = (
    H  = ( 1, 12.00),   # GS98
    He = ( 4, 10.93),   # GS98
    C  = (12,  8.50),   # Caffau+2010
    N  = (14,  7.86),   # Caffau+2009
    O  = (16,  8.76),   # Caffau+2008b  (alpha)
    Ne = (20,  8.08),   # GS98          (alpha)
    Na = (23,  6.33),   # GS98
    Mg = (24,  7.58),   # GS98          (alpha)
    Al = (27,  6.47),   # GS98
    Si = (28,  7.55),   # GS98          (alpha)
    P  = (31,  5.46),   # Caffau+2007b
    S  = (32,  7.16),   # CaffauLudwig+2007  (alpha)
    Ar = (36,  6.40),   # GS98          (alpha)
    Ca = (40,  6.36),   # GS98          (alpha)
    Ti = (48,  5.02),   # GS98          (alpha)
    Fe = (56,  7.52),   # Caffau+2011
    Ni = (58,  6.25),   # GS98
)

const _alpha_elements_parsec = (:O, :Ne, :Mg, :Si, :S, :Ar, :Ca, :Ti)

# Mass contribution of each element: A_i · 10^(ε_i − 12)
_mass(A, ε) = A * exp10(ε - 12.0)

_contributions_parsec = NamedTuple{keys(_bressan2012)}(
    _mass(A, ε) for (A, ε) in values(_bressan2012)
)

Z_total_parsec = sum(v for (el, v) in pairs(_contributions_parsec) if el ∉ (:H, :He))
Z_alpha_parsec = sum(v for (el, v) in pairs(_contributions_parsec) if el ∈ _alpha_elements_parsec)

# f_α = fraction of total metal mass contributed by alpha elements
f_alpha_parsec = Z_alpha_parsec / Z_total_parsec

# The stored value alpha_mass_fraction(PARSECChemistry()) = 0.6723 is f_alpha rounded
# to four decimal places.  We confirm the calculation is consistent with the stored value
# to within half a unit in the last place (atol = 5e-5).
@test round(f_alpha_parsec; digits=4) == alpha_mass_fraction(PARSECChemistry())
