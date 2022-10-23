using BicycleEfficiency
using Documenter


makedocs(
         sitename = "BicycleEfficiency.jl",
         modules  = [BicycleEfficiency],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
    repo="github.com/jotas6/BicycleEfficiency.jl",
)
