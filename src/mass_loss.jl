# Functions to calculate mass loss rates for different types of stars

"""Abstract supertype for all models of stellar mass loss. New models `NewModel <: AbstractMassLoss` should implement `Mdot(m::NewModel, args...)` which should return the mass-loss rate in solar masses per year. A generic method is provided to make subtypes callable -- `(m::AbstractMassLoss)(args...) = Mdot(m, args...)`."""
abstract type AbstractMassLoss{T} end
(m::AbstractMassLoss)(args...) = Mdot(m, args...)

"""
    Bjorklund2021MassLoss(A = -555//100, 
                          B = 79//100, 
                          C = 216//100, 
                          D = -32//100, 
                          E = 6//1, 
                          Zsol = 13//1000) <: AbstractMassLoss
Stellar mass-loss rate model from [Bjorklund2021](@citet), specifically their Equation 20, which reads

```math
\\text{log} \\left( \\dot{M} \\right) = -5.55 + 0.79 \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) + \\left[ 2.16 - 0.32 \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) \\right] \\, \\text{log} \\left( \\frac{L_*}{10^6 \\, L_\\odot} \\right)
```

Note that their model grid assumed ``Z_\\odot = 0.013`` following [Asplund2009](@citet). Internally, assuming the same variable names as used in the call signature above, we compute this as

```math
\\text{log} \\left( \\dot{M} \\right) = A + B \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) + \\left[ C + D \\, \\text{log} \\left( \\frac{ \\text{Z}_* }{ \\text{Z}_\\odot } \\right) \\right] \\, \\left[ \\text{log} \\left( L_* \\right) - E \\right]
```

This model was found to be a better fit to the low-metallicity data of [Telford2024](@citet) and [Hawcroft2024](@citet) than the more commonly used [Vink2001](@citet) model.

Instances are callable with `(Z, logL)` arguments and **return the mass-loss rate in solar masses per year**, where `Z` is metal mass fraction and `logL` is the base-10 logarithm of the star's luminosity in units of solar luminosities.

```jldoctest
julia> using BolometricCorrections: Bjorklund2021MassLoss

julia> model = Bjorklund2021MassLoss();

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

"""
Krticka2025MassLoss(A = -772//100, 
                    B = 149//100, 
                    C = 713//1000, 
                    D = 129//100, 
                    E = 110//100, 
                    T0 = 144//10, 
                    ΔT = 253//100, 
                    Zsol = 13//1000) <: AbstractMassLoss

Stellar mass-loss rate model from [Krticka2025](@citet), specifically their Equation 2, which reads

```math
\\begin{align*}
\\text{log} \\left( \\dot{M} \\right) &= A + B \\, \\text{log} \\left( \\frac{L}{ 10^6 \\, L_\\odot} \\right) + C \\, \\text{log} \\left( \\frac{Z}{Z_\\odot} \\right) + D \\, \\text{log} \\left( \\frac{T}{10^3 \\, K} \\right) + \\\\
& E \\, \\left( \\frac{Z}{Z_\\odot} \\right) \\, \\exp \\left( -\\frac{(T - T_0)^2}{\\Delta T^2} \\right)
\\end{align*}
```

These are line-driven wind models run for OB stars with metallicities down to 0.01 ``Z_\\odot`` with the METUJE code.

Instances are callable with `(Z, logL, Teff)` arguments and **return the mass-loss rate in solar masses per year**, where `Z` is metal mass fraction, `logL` is the base-10 logarithm of the star's luminosity in units of solar luminosities, and `Teff` is the effective temperature of the star in Kelvin.
"""
struct Krticka2025MassLoss{T} <: AbstractMassLoss{T}
    A::T # Additive prefactor
    B::T # Multiplicative prefactor on log(L/10^6 Lsol)
    C::T # Multiplicative prefactor on log(Z/Zsol)
    D::T # Multiplicative prefactor on log(T/1e3 K)
    E::T # Multiplicative prefactor for (Z/Zsol) * exp(-(Teff - T0)^2 / ΔT^2)
    T0::T # kK
    ΔT::T # kK
    Zsol::T # Solar metallicity
end
function Krticka2025MassLoss(A = -772//100, B = 149//100, C = 713//1000, D = 129//100, E = 110//100, T0 = 144//10, ΔT = 253//100, Zsol = 13//1000)
    return Krticka2025MassLoss(promote(A, B, C, D, E, T0, ΔT, Zsol)...)
end
function Mdot(m::Krticka2025MassLoss, Z, logL, Teff)
    # If arguments outside of model grid, return zero mass-loss rate
    if !((0.01 <= Z <= 1) && (10_000 <= Teff <= 45_000))
        return zero(promote_type(typeof(Z), typeof(logL), typeof(Teff)))
    end
    logZ = log10(Z / m.Zsol)
    logMdot = m.A + m.B * (logL - 6) + m.C * logZ + m.D * log10(Teff / 1000) + m.E * (Z / m.Zsol) * exp(-(Teff - m.T0)^2 / m.ΔT^2)
    return exp10(logMdot)
end