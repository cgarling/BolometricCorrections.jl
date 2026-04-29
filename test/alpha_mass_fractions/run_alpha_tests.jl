@testset verbose=true "alpha_mass_fraction derivations" begin
    @safetestset "MISTv1Chemistry (Asplund+2009 protostellar)" include("mist_v1_asplund2009.jl")
    @safetestset "MISTv2Chemistry (Grevesse+1998 photospheric)" include("mist_v2_grevesse1998.jl")
    @safetestset "ATLAS9Chemistry (Grevesse+1998 photospheric)" include("atlas9_grevesse1998.jl")
    @safetestset "PARSECChemistry (Bressan+2012: GS98 + Caffau+2011)" include("parsec_bressan2012.jl")
end
