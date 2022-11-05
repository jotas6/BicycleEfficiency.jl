using BicycleEfficiency
using Test

include("../notebooks/test_data.jl")

N = [50, 21]
ω1 = 50*(2π/60)

@test P1T(Pin, μ, ρ, N, ω1, m, p) > 0 # Test for P1
@test P2T(Pin, N, ω1, μ, rR, γ, m, p) > 0 # Test for P2
@test P3T(Pin, μ, rR, N, ω1, ψ, m, p) > 0 # Test for P3
@test P1T(Pin, μ, ρ, N, ω1, m, p) > P2T(Pin, N, ω1, μ, rR, γ, m, p)
@test P1T(Pin, μ, ρ, N, ω1, m, p) > P3T(Pin, μ, rR, N, ω1, ψ, m, p)
@test PtotalT(Pin, fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, m) > 0 # Test for total power loss
@test P1T(Pin, μ, ρ, N, ω1, m, p) + P2T(Pin, N, ω1, μ, rR, γ, m, p) + P3T(Pin, μ, rR, N, ω1, ψ, m, p) == PtotalT(Pin, fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, m)
@test ηT(Pin, fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, m) > 0.8 # Test for efficiency calculation
@test ηT(Pin, fill(μ, 3), p, ρ, ψ, rR, N, ω1, γ, m) < 1
