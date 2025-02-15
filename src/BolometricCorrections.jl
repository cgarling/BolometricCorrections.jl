module BolometricCorrections

import CSV
using TypedTables: Table, columnnames

# exports from top-level module
# const test = "asdf" # This is exported at the module level
# export test

include("YBC/YBC.jl")
using .YBC
# exports from YBC
include("MIST/MIST.jl")
using .MIST
# exports from MIST

end # module
