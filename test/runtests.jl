import BolometricCorrections as BC
using SafeTestsets
using Test

# Run doctests first
using Documenter: DocMeta, doctest
DocMeta.setdocmeta!(BC, :DocTestSetup, :(using BolometricCorrections); recursive=true)
doctest(BC)

# Run tests for MIST submodule
include(joinpath("MIST", "run_mist_tests.jl"))
include(joinpath("YBC", "run_ybc_tests.jl"))
