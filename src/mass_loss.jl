# Functions to calculate mass loss rates for different types of stars

"""Abstract supertype for all models of stellar mass loss."""
abstract type AbstractMassLoss{T} end

"""
    Bjorklund2021MassLoss(A = -555//100, B = 79//100, C = 216//100, D = -32//100, E = 6//1, Zsol = 13//1000)
Stellar mass-loss rate model from [Bjorklund2021](@citet), specifically their Equation 20, which reads

```math
\\text{log} \\left( \\dot{M} \\right) = -5.55 + 0.79 \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) + \\left[ 2.16 - 0.32 \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) \\right] \\, \\text{log} \\left( \\frac{L_*}{10^6 \\, L_\\odot} \\right)
```

Note that their model grid assumed ``Z_\\odot = 0.013`` following [Asplund2009](@citet). Internally, assuming the same variable names as used in the call signature above, we compute this as

```math
\\text{log} \\left( \\dot{M} \\right) = A + B \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) + \\left[ C + D \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) \\right] \\, \\left[ \\text{log} \\left( L_* \\right) - E \\right]
```

This model was found to be a better fit to the low-metallicity data of [Telford2024](@citet) than the more commonly used [Vink2001](@citet) model.

Instances are callable with `(Z, logL)` arguments and **return the mass-loss rate in solar masses per year**, where `Z` is metal mass fraction and `logL` is the base-10 logarithm of the star's luminosity in units of solar luminosities.

```jldoctest
julia> using BolometricCorrections: Bjorklund2021MassLoss

julia> model = Bjorklund2021MassLoss()

julia> model(1e-3, 5) ≈ 1.1311569779109792e-9 # Mdot in solar mass per year
true
```
"""
struct Bjorklund2021MassLoss{T} <: AbstractMassLoss{T}
    A::T # Additive prefactor, -5.55
    B::T # Multiplicative prefactor on log(Z)
    C::T
    D::T
    E::T # luminosity scaling is log10(L_* / L_ref) -- this value is log10(L_ref)
    Zsol::T
end
function Bjorklund2021MassLoss(A = -555//100, B = 79//100, C = 216//100, D = -32//100, E = 6//1, Zsol = 13//1000)
    return Bjorklund2021MassLoss(promote(A, B, C, D, E, Zsol)...)
end
function Mdot(m::Bjorklund2021MassLoss, Z, logL)
    logZ = log10(Z / m.Zsol)
    logMdot = m.A + m.B * logZ + (m.C + m.D * logZ) * (logL - m.E)
    return exp10(logMdot)
end
(m::Bjorklund2021MassLoss)(Z, logL) = Mdot(m, Z, logL)