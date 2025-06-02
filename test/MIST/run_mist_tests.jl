@testset verbose=true "MIST Submodule" begin
    @safetestset "select_subtable" include("select_subtable.jl")
    @safetestset "chemistry" include("mist_chemistry.jl")
    @safetestset "zeropoints" include("mist_zpt_test.jl")
end
