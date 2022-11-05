using BicycleEfficiency
using Test

include("../notebooks/test_data.jl")

N = [50, 21]
ω1 = 50*(2π/60)

@test P1(μ, ρ, N, ω1, T0) > 0 # Test for P1
@test P2(N, ω1, μ, rR, γ, T0) > 0 # Test for P2
@test P3(μ, rR, N, ω1, ψ, T0) > 0 # Test for P3
@test P1(μ, ρ, N, ω1, T0) > P2(N, ω1, μ, rR, γ, T0)
@test P1(μ, ρ, N, ω1, T0) > P3(μ, rR, N, ω1, ψ, T0)
@test Ptotal(fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, T0) > 0 # Test for total power loss
@test P1(μ, ρ, N, ω1, T0) + P2(N, ω1, μ, rR, γ, T0) + P3(μ, rR, N, ω1, ψ, T0) == Ptotal(fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, T0)
@test η(fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, T0) > 0.8 # Test for efficiency calculation
@test η(fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, T0) < 1
