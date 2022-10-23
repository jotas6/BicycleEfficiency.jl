"""
    P1(μ1::Number, ρ::Number, T0::Number, N::Vector{Int}, ω::Number)::Float64

Compute the power loss due to friction between pins and bushings.

# Arguments
- `μ1::Number`: friction coefficient between the pin and the bushing.
- `ρ::Number`: bushing radius.
- `T0::Number`: free chain tension.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
"""
function P1(μ1::Number, ρ::Number, T0::Number, N::Vector{Int}, ω::Number)::Float64

    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)

    p1 = μ1*ρ*(π/2)*T0

    num = @. 1 + tan(α/2)/tan(ϕ/2)
    den = @. 1 - tan(α/2)/tan(ϕ/2)

    p2 = @. sin(ϕ)log(abs(num/den))

    Wf = p1*sum(p2)

    Pf = N[1]*ω*Wf/2π

    return Pf

end

"""
    P2(N::Vector{Int}, ω::Number, μ2::Number, T0::Number, r0::Number, γ::Number)::Float64

Compute the power loss due to chain offset

# Arguments
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence
- `μ2::Number`: friction coefficient between the chain and the sprocket teeth
- `T0::Number`: free chain tension
- `r0::Number`: radius of contact during offset operation
- `γ::Number`: chain offset angle
"""
function P2(N::Vector{Int}, ω::Number, μ2::Number, T0::Number, r0::Number, γ::Number)::Float64

    p1 = N[1]*ω*μ2
    p2 = T0*r0*sin(γ)
    p3 = sum(1 ./ N)

    Pf = p1*p2*p3

    return Pf

end

"""
    P3(μ3::Number, T0::Number, rR::Number, N::Vector{Int}, ω::Number, ψ::Number)::Float64

Compute the power loss due to interaction between rollers and sprocket teeth

# Arguments
- `μ3::Number`: friction coefficient between roller and sprocket teeth.
- `T0::Number`: free chain tension.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence
- `ψ::Number`: absolute roller rotation angle.
"""
function P3(μ3::Number, T0::Number, rR::Number, N::Vector{Int}, ω::Number, ψ::Number)::Float64

    δ = ψ*rR

    ϕ = @. (30 - 120/N)*(π/180)

    p1 = μ3*T0*rR*N[1]*ω/(4π)

    p21 = ((N[1]*ψ/(π - 1))*δ + (N[2]*ψ/(π + 1)*δ))

    p221 = @. ((2π)/N)cos(ϕ)
    p222 = @. sin(ϕ)log(cos((2π)/N) + sin((2π)/N)*cot(ϕ))

    p22 = p221 - p222

    p2 = @. p21*p22

    Pf = p1*sum(p2)

    return Pf

end


 """
    Ptotal(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,
           T0::Number, N::Vector{Int}, ω::Number, γ::Number)::Float64

Compute the power loss due to all 3 sources (friction between pins and bushings,
chain misalignment, interaction between rollers and sprocket teeth)

It is calculated as the sum of the losses from all sources, which are implemented
separately in P1, P2 and P3

# Arguments
- `μ::Vector{Float64}`: vector containing the friction coefficients for cases 1, 2 and 3 (in that order).
- `p::Number`: chain pitch.
- `ρ::Number`: bushing radius.
- `rR::Number`: roller radius.
- `T0::Number`: free chain tension.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `γ::Number`: chain offset angle
"""
function Ptotal(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,
                T0::Number, N::Vector{Int}, ω::Number, γ::Number)::Float64

    δ = ψ*rR
    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)

    P1_value = P1(μ[1], ρ, T0, N, ω)
	P2_value = P2(N, ω, μ[2], T0, rR, γ)
	P3_value = P3(μ[3], T0, rR, N, ω, ψ)

	return P1_value + P2_value + P3_value

end


"""
    η(μ::Vector{Float64}, p::Float64, ρ::Float64, ψ::Float64, rR::Float64,
      T0::Float64, N::Vector{Int}, ω::Float64)::Float64

Compute the power transmission efficiency considering only frictional losses

# Arguments
- `μ::Vector{Float64}`: friction coefficients for cases 1, 2 and 3 (in that order).
- `p::Number`: chain pitch.
- `ρ::Number`: bushing radius.
- `ψ::Number`: absolute roller rotation angle.
- `rR::Number`: roller radius.
- `T0::Number`: free chain tension.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
"""
function η(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,
           T0::Number, N::Vector{Int}, ω::Number, γ::Number)::Float64

    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)
    δ = ψ*rR

	Ptotal_value = Ptotal(μ, p, ρ, ψ, rR, T0, N, ω, γ)
	Pin = p*T0*ω/α[1]

	return 1 - Ptotal_value/Pin

end
