using BicycleEfficiency
using Test

include("../notebooks/test_data.jl")

N = [50, 21]
ω1 = 50*(2π/60)

@testset "BicycleEfficiency.jl" begin
    @test P1(μ, ρ, T0, N, ω1) > 0 # Test for P1
    @test P2(N, ω1, μ, T0, rR, γ) > 0 # Test for P2
    @test P3(μ, T0, rR, N, ω1, ψ) > 0 # Test for P3
    @test P1(μ, ρ, T0, N, ω1) > P2(N, ω1, μ, T0, rR, γ)
    @test P1(μ, ρ, T0, N, ω1) > P3(μ, T0, rR, N, ω1, ψ)
    @test Ptotal(fill(μ, 3), p, ρ, ψ, rR, T0, N, ω1, γ) > 0 # Test for total power loss
    @test P1(μ, ρ, T0, N, ω1) + P2(N, ω1, μ, T0, rR, γ) + P3(μ, T0, rR, N, ω1, ψ) == Ptotal(fill(μ, 3), p, ρ, ψ, rR, T0, N, ω1, γ)
    @test η(fill(μ, 3), p, ρ, ψ, rR, T0, N, ω1, γ) > 0.8 # Test for efficiency calculation
    @test η(fill(μ, 3), p, ρ, ψ, rR, T0, N, ω1, γ) < 1
end
