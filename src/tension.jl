"""
    Tf(Pin::Number, ω1::Number, p::Number, N1::Int)::Float64

Compute the chain tension due to the force applied by the sprocket

# Arguments
- `Pin::Number`: power input
- `ω::Number`: pedaling cadence
- `p::Number`: chain pitch
- `N1::Int`: number of teeth on the front sprocket

# Examples
```julia-repl
julia> Tf(150, 80*(2π/60), 12.22e-3, 48)
191.79623567921442
"""
function Tf(Pin::Number, ω::Number, p::Number, N1::Int)::Float64

    return (2π*Pin)/(ω*p*N1)

end

"""
    Tc(m::Number, ω::Number, p::Number, N1::Int)::Float64

Compute the chain tension due to centrifugal force

# Arguments
- `m::Number`: chain linear density (mass/length)
- `ω::Number`: pedaling cadence
- `p::Number`: chain pitch
- `N1::Int`: number of teeth on the front sprocket

# Examples
```julia-repl
julia> Tc(1.52, 80*(2π/60), 12.22e-3, 48)
0.929706672128
```
"""
function Tc(m::Number, ω::Number, p::Number, N1::Int)::Float64

    return (m*(ω*p*N1)^2)/(4π^2)

end

"""
    T(Pin::Number, ω::Number, p::Number, N1::Int, m::Number)::Float64

Compute total chain tension, which is the sum of the tension caused
by the force applied by the sprocket and the tension due to centrifugal
forces

# Arguments
- `Pin::Number`: power input
- `ω::Number`: pedaling cadence
- `p::Number`: chain pitch
- `N1::Int`: number of teeth on the front sprocket
- `m::Number`: chain linear density (mass/length)

# Examples
```julia-repl
julia> T(150, 80*(2π/60), 12.22e-3, 48, 1.52)
192.7259423513424
```
"""
function T(Pin::Number, ω::Number, p::Number, N1::Int, m::Number)::Float64

    Tf_value = Tf(Pin, ω, p, N1)
    Tc_value = Tc(m, ω, p, N1)

    return Tf_value + Tc_value
end
