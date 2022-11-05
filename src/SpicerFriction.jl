
"""
    P1(μ1::Number, ρ::Number, N::Vector{Int}, ω::Number, T0::Number)::Float64

Compute the power loss due to friction between pins and bushings.

# Arguments
- `μ1::Number`: friction coefficient between the pin and the bushing.
- `ρ::Number`: bushing radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `T0::Number`: free chain tension

# Examples
```julia-repl
julia> P1(0.09, 0.00175968, 300, [48, 24], 80*(2π/60), 200)
1.8756733110309638
```
"""
function P1(μ1::Number, ρ::Number, N::Vector{Int}, ω::Number, T0::Number)::Float64

    return N[1]*ω*μ1*ρ*(π/2)*T0*sum(1 ./ N)

end

"""
    P2(N::Vector{Int}, ω::Number, μ2::Number, r0::Number, γ::Number, T0::Number)::Float64

Compute the power loss due to chain offset.

# Arguments
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence
- `μ2::Number`: friction coefficient between the chain and the sprocket teeth
- `r0::Number`: radius of contact during offset operation
- `γ::Number`: chain offset angle
- `T0::Number`: free chain tension

# Examples
```julia-repl
julia> P2(150, [48, 24], 80*(2pi/60), 0.09, 300, 0.00249288, π/180, 1.52)
0.02952298838057144
```
"""
function P2(N::Vector{Int}, ω::Number, μ2::Number, r0::Number, γ::Number, T0::Number)::Float64

    p1 = N[1]*ω*μ2
    p2 = T0*r0*sin(γ)
    p3 = sum(1 ./ N)

    Pf = p1*p2*p3

    return Pf

end

"""
    P3(μ3::Number, rR::Number, N::Vector{Int}, ω::Number, ψ::Number, T0::Number)::Float64

Compute the power loss due to interaction between rollers and sprocket teeth.

# Arguments
- `μ3::Number`: friction coefficient between roller and sprocket teeth.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence
- `ψ::Number`: absolute roller rotation angle.
- `T0::Number`: free chain tension

# Examples
```julia-repl
julia> P3(0.09, 300, 0.00249288, [48, 24], 80*(2pi/60), π/2)
0.02840827017479334
```
"""
function P3(μ3::Number, rR::Number, N::Vector{Int}, ω::Number, ψ::Number, T0::Number)::Float64

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
           N::Vector{Int}, ω::Number, γ::Number, T0::Number)::Float64

Compute the power loss due to all 3 sources (friction between pins and bushings,
chain misalignment and interaction between rollers and sprocket teeth).

It is calculated as the sum of the losses from all sources, which are implemented
separately in P1, P2 and P3

# Arguments
- `μ::Vector{Float64}`: vector containing the friction coefficients for cases 1, 2 and 3 (in that order).
- `p::Number`: chain pitch.
- `ρ::Number`: bushing radius.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `γ::Number`: chain offset angle
- `T0::Number`: free chain tension

# Examples
```julia-repl
julia> Ptotal(fill(0.09, 3), 0.01222, 0.00175968, π/2, 0.00249288, 300, [48, 24], 80*(2π/60), π/180)
1.9336045695863286
```
"""
function Ptotal(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,
                N::Vector{Int}, ω::Number, γ::Number, T0::Number)::Float64

    δ = ψ*rR
    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)

    P1_value = P1(μ[1], ρ, N, ω, T0)
	P2_value = P2(N, ω, μ[2], rR, γ, T0)
	P3_value = P3(μ[3], rR, N, ω, ψ, T0)

	return P1_value + P2_value + P3_value

end


"""
    η(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number,
      rR::Number, N::Vector{Int}, ω::Number, γ::Number, T0::Number)::Float64

Compute the power transmission efficiency considering only frictional losses.

# Arguments
- `μ::Vector{Float64}`: friction coefficients for cases 1, 2 and 3 (in that order).
- `p::Number`: chain pitch.
- `ρ::Number`: bushing radius.
- `ψ::Number`: absolute roller rotation angle.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `γ::Number`: chain offset angle
- `T0::Number`: free chain tension

# Examples
```julia-repl
julia> η(fill(0.09, 3), 0.01222, 0.00175968, π/2, 0.00249288, 300, [48, 24], 80*(2π/60), π/180, 200)
0.9917587093835826
```
"""
function η(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number,
           rR::Number, N::Vector{Int}, ω::Number, γ::Number, T0::Number)::Float64

    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)
    δ = ψ*rR

	Ptotal_value = Ptotal(μ, p, ρ, ψ, rR, N, ω, γ, T0)
	Pin = p*T0*ω/α[1]

	return 1 - Ptotal_value/Pin

end

################################################################################
# Functions with chain tension model

include("tension.jl")

"""
    P1T(Pin::Number, μ1::Number, ρ::Number, N::Vector{Int}, ω::Number,
       m::Number, p::Number)::Float64

Compute the power loss due to friction between pins and bushings, using chain tension
model.

# Arguments
- `Pin::Number`: power input
- `μ1::Number`: friction coefficient between the pin and the bushing.
- `ρ::Number`: bushing radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `m::Number`: chain linear density (mass/length)
- `p::Number`: chain pitch

# Examples
```julia-repl
julia> P1(150, 0.09, 0.00175968, 300, [48, 24], 80*(2π/60), 1.52)
1.8756733110309638
```
"""
function P1T(Pin::Number, μ1::Number, ρ::Number, N::Vector{Int}, ω::Number, m::Number, p::Number)::Float64

    T0 = T(Pin, ω, p, N[1], m)

    return N[1]*ω*μ1*ρ*(π/2)*T0*sum(1 ./ N)

end

"""
    P2T(Pin::Number, N::Vector{Int}, ω::Number, μ2::Number, r0::Number,
       γ::Number, m::Number)::Float64

Compute the power loss due to chain offset, using chain tension model.

# Arguments
- `Pin::Number`: power input
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence
- `μ2::Number`: friction coefficient between the chain and the sprocket teeth
- `r0::Number`: radius of contact during offset operation
- `γ::Number`: chain offset angle
- `m::Number`: chain linear density (mass/length)
- `p::Number`: chain pitch

# Examples
```julia-repl
julia> P2(150, [48, 24], 80*(2pi/60), 0.09, 300, 0.00249288, π/180, 1.52)
0.02952298838057144
```
"""
function P2T(Pin::Number, N::Vector{Int}, ω::Number, μ2::Number, r0::Number,
            γ::Number, m::Number, p::Number)::Float64

    T0 = T(Pin, ω, p, N[1], m)

    p1 = N[1]*ω*μ2
    p2 = T0*r0*sin(γ)
    p3 = sum(1 ./ N)

    Pf = p1*p2*p3

    return Pf

end

"""
    P3T(μ3::Number, T0::Number, rR::Number, N::Vector{Int}, ω::Number, ψ::Number)::Float64

Compute the power loss due to interaction between rollers and sprocket teeth, using chain tension
model

# Arguments
- `Pin::Number`: power input
- `μ3::Number`: friction coefficient between roller and sprocket teeth.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence
- `ψ::Number`: absolute roller rotation angle.
- `m::Number`: chain linear density (mass/length)
- `p::Number`: chain pitch

# Examples
```julia-repl
julia> P3(0.09, 300, 0.00249288, [48, 24], 80*(2pi/60), π/2)
0.02840827017479334
```
"""
function P3T(Pin::Number, μ3::Number, rR::Number, N::Vector{Int}, ω::Number,
            ψ::Number, m::Number, p::Number)::Float64

    T0 = T(Pin, ω, p, N[1], m)

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
    PtotalT(Pin::Number, μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,
           N::Vector{Int}, ω::Number, γ::Number, m::Number)::Float64

Compute the power loss due to all 3 sources (friction between pins and bushings,
chain misalignment and interaction between rollers and sprocket teeth), using chain
tension model.

It is calculated as the sum of the losses from all sources, which are implemented
separately in P1, P2 and P3

# Arguments
- `Pin::Number`: power input
- `μ::Vector{Float64}`: vector containing the friction coefficients for cases 1, 2 and 3 (in that order).
- `p::Number`: chain pitch.
- `ρ::Number`: bushing radius.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `γ::Number`: chain offset angle
- `m::Number`: chain linear density (mass/length)

# Examples
```julia-repl
julia> Ptotal(fill(0.09, 3), 0.01222, 0.00175968, π/2, 0.00249288, 300, [48, 24], 80*(2π/60), π/180)
1.9336045695863286
```
"""
function PtotalT(Pin::Number, μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,
                N::Vector{Int}, ω::Number, γ::Number, m::Number)::Float64

    T0 = T(Pin, ω, p, N[1], m)

    δ = ψ*rR
    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)

    P1_value = P1T(Pin, μ[1], ρ, N, ω, m, p)
	P2_value = P2T(Pin, N, ω, μ[2], rR, γ, m, p)
	P3_value = P3T(Pin, μ[3], rR, N, ω, ψ, m, p)

	return P1_value + P2_value + P3_value

end


"""
    ηT(Pin::Number, rμ::Vector{Float64}, p::Number, ρ::Number, ψ::Number,
      rR::Number, N::Vector{Int}, ω::Number, γ::Number, m::Number)::Float64

Compute the power transmission efficiency considering only frictional losses,
using chain tension model

# Arguments
- `Pin::Number`: power input
- `μ::Vector{Float64}`: friction coefficients for cases 1, 2 and 3 (in that order).
- `p::Number`: chain pitch.
- `ρ::Number`: bushing radius.
- `ψ::Number`: absolute roller rotation angle.
- `rR::Number`: roller radius.
- `N::Vector{Int}`: number of teeth on the front and rear sprockets (in that order).
- `ω::Number`: pedaling cadence.
- `γ::Number`: chain offset angle
- `m::Number`: chain linear density (mass/length)

# Examples
```julia-repl
julia> η(150, fill(0.09, 3), 0.01222, 0.00175968, π/2, 0.00249288, 300, [48, 24], 80*(2π/60), π/180, 1.52)
0.9917587093835826
```
"""
function ηT(Pin::Number, μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number,
           rR::Number, N::Vector{Int}, ω::Number, γ::Number, m::Number)::Float64

    T0 = T(Pin, ω, p, N[1], m)

    α = @. (360/N)*(π/180)
    ϕ = @. (30 - 120/N)*(π/180)
    δ = ψ*rR

	Ptotal_value = PtotalT(Pin, μ, p, ρ, ψ, rR, N, ω, γ, m)

	return 1 - Ptotal_value/Pin

end
