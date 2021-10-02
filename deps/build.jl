append!(empty!(LOAD_PATH), Base.DEFAULT_LOAD_PATH)
using Pkg
exampledir = joinpath(@__DIR__, "..", "examples")
Pkg.activate(exampledir)
Pkg.add("Literate")
include(joinpath(exampledir, "generate_notebooks.jl"))
