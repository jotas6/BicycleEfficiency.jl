using BicycleEfficiency
using Test

include("../notebooks/test_data.jl")

N = [50, 21]
ω1 = 50*(2π/60)

@testset "BicycleEfficiency.jl" begin
    begin include("SpicerFrictionTest.jl") end
    begin include("SpicerFrictionTestTension.jl") end
end
