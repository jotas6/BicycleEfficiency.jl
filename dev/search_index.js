var documenterSearchIndex = {"docs":
[{"location":"#BicycleEfficiency.jl","page":"Home","title":"BicycleEfficiency.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package to calculate the efficiency of a bicycle's chain drive","category":"page"},{"location":"#Frictional-losses","page":"Home","title":"Frictional losses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"P1\nP2\nP3\nPtotal\nη","category":"page"},{"location":"#BicycleEfficiency.P1","page":"Home","title":"BicycleEfficiency.P1","text":"P1(μ1::Number, ρ::Number, T0::Number, N::Vector{Int}, ω::Number)::Float64\n\nCompute the power loss due to friction between pins and bushings.\n\nArguments\n\nμ1::Number: friction coefficient between the pin and the bushing.\nρ::Number: bushing radius.\nT0::Number: free chain tension.\nN::Vector{Int}: number of teeth on the front and rear sprockets (in that order).\nω::Number: pedaling cadence.\n\nExamples\n\n\njulia> P1(0.09, 0.00175968, 300, [48, 24], 80*(2π/60))\n3.970776399955861\n\n\n\n\n\n\n","category":"function"},{"location":"#BicycleEfficiency.P2","page":"Home","title":"BicycleEfficiency.P2","text":"P2(N::Vector{Int}, ω::Number, μ2::Number, T0::Number, r0::Number, γ::Number)::Float64\n\nCompute the power loss due to chain offset\n\nArguments\n\nN::Vector{Int}: number of teeth on the front and rear sprockets (in that order).\nω::Number: pedaling cadence\nμ2::Number: friction coefficient between the chain and the sprocket teeth\nT0::Number: free chain tension\nr0::Number: radius of contact during offset operation\nγ::Number: chain offset angle\n\nExamples\n\n\njulia> P2([48, 24], 80*(2pi/60), 0.09, 300, 0.00249288, π/180)\n0.02952298838057144\n\n\n\n\n\n\n","category":"function"},{"location":"#BicycleEfficiency.P3","page":"Home","title":"BicycleEfficiency.P3","text":"P3(μ3::Number, T0::Number, rR::Number, N::Vector{Int}, ω::Number, ψ::Number)::Float64\n\nCompute the power loss due to interaction between rollers and sprocket teeth\n\nArguments\n\nμ3::Number: friction coefficient between roller and sprocket teeth.\nT0::Number: free chain tension.\nrR::Number: roller radius.\nN::Vector{Int}: number of teeth on the front and rear sprockets (in that order).\nω::Number: pedaling cadence\nψ::Number: absolute roller rotation angle.\n\nExamples\n\n\njulia> P3(0.09, 300, 0.00249288, [48, 24], 80*(2pi/60), π/2)\n0.02840827017479334\n\n\n\n\n\n\n","category":"function"},{"location":"#BicycleEfficiency.Ptotal","page":"Home","title":"BicycleEfficiency.Ptotal","text":"Ptotal(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,\n       T0::Number, N::Vector{Int}, ω::Number, γ::Number)::Float64\n\nCompute the power loss due to all 3 sources (friction between pins and bushings, chain misalignment, interaction between rollers and sprocket teeth)\n\nIt is calculated as the sum of the losses from all sources, which are implemented separately in P1, P2 and P3\n\nArguments\n\nμ::Vector{Float64}: vector containing the friction coefficients for cases 1, 2 and 3 (in that order).\np::Number: chain pitch.\nρ::Number: bushing radius.\nrR::Number: roller radius.\nT0::Number: free chain tension.\nN::Vector{Int}: number of teeth on the front and rear sprockets (in that order).\nω::Number: pedaling cadence.\nγ::Number: chain offset angle\n\nExamples\n\n\njulia> Ptotal(fill(0.09, 3), 0.01222, 0.00175968, π/2, 0.00249288, 300, [48, 24], 80*(2π/60), π/180)\n4.028707658511226\n\n\n\n\n\n\n","category":"function"},{"location":"#BicycleEfficiency.η","page":"Home","title":"BicycleEfficiency.η","text":"η(μ::Vector{Float64}, p::Number, ρ::Number, ψ::Number, rR::Number,\n  T0::Number, N::Vector{Int}, ω::Number, γ::Number)::Float64\n\nCompute the power transmission efficiency considering only frictional losses\n\nArguments\n\nμ::Vector{Float64}: friction coefficients for cases 1, 2 and 3 (in that order).\np::Number: chain pitch.\nρ::Number: bushing radius.\nψ::Number: absolute roller rotation angle.\nrR::Number: roller radius.\nT0::Number: free chain tension.\nN::Vector{Int}: number of teeth on the front and rear sprockets (in that order).\nω::Number: pedaling cadence.\nγ::Number: chain offset angle\n\nExamples\n\n\njulia> η(fill(0.09, 3), 0.01222, 0.00175968, π/2, 0.00249288, 300, [48, 24], 80*(2π/60), π/180)\n0.9828290896987895\n\n\n\n\n\n\n","category":"function"}]
}
