push!(LOAD_PATH,"../src/")

using Documenter, QuasiGeostrophy

makedocs(
    modules = [QuasiGeostrophy],
    pages = [
        "Home" => "index.md",
    ],
    sitename = "QuasiGeostrophy"
)

deploydocs(repo = "github.com/sandreza/QuasiGeostrophy.git")

#=
using Documenter, QuasiGeostrophy, Literate


If Literate is included then the following must be added to the .travis.yml file
julia --project=docs/ -e 'using Pkg; Pkg.instantiate();
                    Pkg.add("Documenter");
                    Pkg.add("Literate");
                    Pkg.add("Plots");
                    Pkg.develop(PackageSpec(path=pwd()))'
specifically the lines Pkg.add("Documenter"); Pkg.add("Literate"); Pkg.add("Plots");
=#
# doesn't necessarily need to be generated on the fly
#=
const examples_directory = pwd() * "/docs/example/"
const output_directory = pwd() * "/docs/src/generated/"
# read the examples in the directory
examples = readdir(examples_directory)

for example in examples
    example_filepath = examples_directory * example
    Literate.markdown(example_filepath, output_directory, documenter = true)
end



makedocs(
    modules = [QuasiGeostrophy],
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "first example" => "generated/literate_example.md"
        ],
        "Function Index" => "function_index.md"
    ],
    sitename = "QuasiGeostrophy"
)

deploydocs(repo = "github.com/sandreza/QuasiGeostrophy.git")
=#
