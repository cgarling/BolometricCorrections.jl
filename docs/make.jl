using Documenter
using BolometricCorrections
import BolometricCorrections

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

# We check link validity with `linkcheck=true`, but we don't want this to fail the build
# so we add `:linkcheck` to `warnonly`. Additionally, we are setting
# `modules = [Bolometriccorrections]` so a warning will be raised if any inline
# documentation strings are not included in the document. In v1.0 of Documenter, this warning
# will raise an error and prevent running. By adding `:missing_docs` to `warnonly`, we will
# see these warnings but they will not raise an error.

# DocMeta.setdocmeta!(BolometricCorrections, :DocTestSetup, :(using BolometricCorrections); recursive=true)

makedocs(
    sitename = "BolometricCorrections.jl",
    modules = [BolometricCorrections],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                             size_threshold_warn = 409600, # v1.0.0 default: 102400 (bytes)
                             size_threshold = 819200,      # v1.0.0 default: 204800 (bytes)
                             example_size_threshold=0),    # Write all @example to file
    authors = "Chris Garling",
    pages = ["index.md",
             "MIST.md",
                  # "Internals" => ["fitting/internals.md",
             #                 "fitting/kernels.md"]],
             "api.md",
             "doc_index.md"],
    doctest = false,
    linkcheck = true,
    warnonly = [:missing_docs, :linkcheck]
)

deploydocs(;
    repo = "github.com/cgarling/BolometricCorrections.jl.git",
    versions = ["stable" => "v^", "v#.#"],
    push_preview=true,
)
