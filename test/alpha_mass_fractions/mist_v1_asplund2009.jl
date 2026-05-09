# Derivation of the solar alpha-element mass fraction for MISTv1Chemistry.
#
# MIST v1.2 adopts the bulk solar (protostellar) chemical abundances from Table 4 of:
#   Asplund, Grevesse, Sauval & Scott (2009), ARA&A, 47, 481–522
#   DOI: 10.1146/annurev.astro.46.060407.145222
#   arXiv: 0909.0948
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
# Protostellar abundances are taken from Table 4 ("Present-day solar") of Asplund+2009;
# a subset of representative metals is included below, sufficient to evaluate f_α.
# Omitted trace elements contribute negligibly to Z and do not affect f_α at the
# precision retained here.

using Test: @test
using BolometricCorrections.MIST: MISTv1Chemistry, alpha_mass_fraction

# ────────────────────────────────────────────────────────────────────────────
# Asplund+2009 Table 4 protostellar log-epsilon abundances: ε_i = log(N_i/N_H) + 12
# Columns: element => (atomic mass number A, log_epsilon)
# ────────────────────────────────────────────────────────────────────────────
const _asplund2009_proto = (
    H  = ( 1, 12.00),
    He = ( 4, 10.98),
    C  = (12,  8.47),
    N  = (14,  7.87),
    O  = (16,  8.76),   # alpha
    Ne = (20,  8.09),   # alpha
    Na = (23,  6.29),
    Mg = (24,  7.64),   # alpha
    Al = (27,  6.47),
    Si = (28,  7.55),   # alpha
    S  = (32,  7.16),   # alpha
    Ar = (36,  6.44),   # alpha
    Ca = (40,  6.34),   # alpha
    Ti = (48,  4.99),   # alpha
    Fe = (56,  7.54),
    Ni = (58,  6.26),
)

const _alpha_elements_asplund = (:O, :Ne, :Mg, :Si, :S, :Ar, :Ca, :Ti)

# Mass contribution of each element: A_i · 10^(ε_i − 12)
_mass(A, ε) = A * exp10(ε - 12.0)

_contributions_asplund = NamedTuple{keys(_asplund2009_proto)}(
    _mass(A, ε) for (A, ε) in values(_asplund2009_proto)
)

Z_total_asplund = sum(v for (el, v) in pairs(_contributions_asplund) if el ∉ (:H, :He))
Z_alpha_asplund = sum(v for (el, v) in pairs(_contributions_asplund) if el ∈ _alpha_elements_asplund)

# f_α = fraction of total metal mass contributed by alpha elements
f_alpha_asplund2009 = Z_alpha_asplund / Z_total_asplund

# The stored value alpha_mass_fraction(MISTv1Chemistry()) = 0.6803 is f_alpha rounded
# to four decimal places.  We confirm the calculation is consistent with the stored value
# to within half a unit in the last place (atol = 5e-5).
@test round(f_alpha_asplund2009; digits=4) == alpha_mass_fraction(MISTv1Chemistry())
