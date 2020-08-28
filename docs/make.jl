push!(LOAD_PATH,"../src/")

using Documenter, QuasiGeostrophy


equations = ["Home" => "equations.md",
             "2-Layer QG" => "2_layers.md",
             "N-Layer QG" => "n_layers.md",
             ]

makedocs(
    modules = [QuasiGeostrophy],
    pages = [
        "Home" => "index.md",
        "Equations" => equations,
    ],
    sitename = "QuasiGeostrophy",
    format = Documenter.HTML(collapselevel = 1),
)

deploydocs(repo = "github.com/sandreza/QuasiGeostrophy.git")
