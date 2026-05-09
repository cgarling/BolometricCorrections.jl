# Derivation of the solar alpha-element mass fraction for MISTv2Chemistry.
#
# MIST v2.5 adopts the solar photospheric abundances from Table 1 of:
#   Grevesse, N. & Sauval, A. J. (1998), Space Science Reviews, 85, 161–174
#   DOI: 10.1023/A:1005161325181
#   ADS: https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G
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
# Photospheric abundances are taken from Table 1 of Grevesse & Sauval (1998); a subset of
# representative metals is included below, sufficient to evaluate f_α.
# Omitted trace elements contribute negligibly to Z and do not affect f_α at the
# precision retained here.

using Test: @test
using BolometricCorrections.MIST: MISTv2Chemistry, alpha_mass_fraction

# ────────────────────────────────────────────────────────────────────────────
# Grevesse & Sauval (1998) Table 1 photospheric log-epsilon abundances: ε_i = log(N_i/N_H) + 12
# Columns: element => (atomic mass number A, log_epsilon)
# ────────────────────────────────────────────────────────────────────────────
const _gs98 = (
    H  = ( 1, 12.00),
    He = ( 4, 10.93),
    C  = (12,  8.52),
    N  = (14,  7.92),
    O  = (16,  8.83),   # alpha
    Ne = (20,  8.08),   # alpha
    Na = (23,  6.33),
    Mg = (24,  7.58),   # alpha
    Al = (27,  6.47),
    Si = (28,  7.55),   # alpha
    P  = (31,  5.45),
    S  = (32,  7.33),   # alpha
    Ar = (36,  6.40),   # alpha
    Ca = (40,  6.36),   # alpha
    Ti = (48,  5.02),   # alpha
    Fe = (56,  7.50),
    Ni = (58,  6.25),
)

const _alpha_elements_gs98 = (:O, :Ne, :Mg, :Si, :S, :Ar, :Ca, :Ti)

# Mass contribution of each element: A_i · 10^(ε_i − 12)
_mass(A, ε) = A * exp10(ε - 12.0)

_contributions_gs98 = NamedTuple{keys(_gs98)}(
    _mass(A, ε) for (A, ε) in values(_gs98)
)

Z_total_gs98 = sum(v for (el, v) in pairs(_contributions_gs98) if el ∉ (:H, :He))
Z_alpha_gs98 = sum(v for (el, v) in pairs(_contributions_gs98) if el ∈ _alpha_elements_gs98)

# f_α = fraction of total metal mass contributed by alpha elements
f_alpha_gs98 = Z_alpha_gs98 / Z_total_gs98

# The stored value alpha_mass_fraction(MISTv2Chemistry()) = 0.6911 is f_alpha rounded
# to four decimal places.  We confirm the calculation is consistent with the stored value
# to within half a unit in the last place (atol = 5e-5).
@test round(f_alpha_gs98; digits=4) == alpha_mass_fraction(MISTv2Chemistry())
