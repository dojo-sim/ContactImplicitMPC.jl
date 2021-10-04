append!(empty!(LOAD_PATH), Base.DEFAULT_LOAD_PATH)
using Pkg

################################################################################
# Generate notebooks
################################################################################
exampledir = joinpath(@__DIR__, "..", "examples")
Pkg.activate(exampledir)
Pkg.add("Literate")
include(joinpath(exampledir, "generate_notebooks.jl"))

# ################################################################################
# # Code generation
# ################################################################################
# # Add minimal set of packages
# Pkg.add(PackageSpec(name = "FileIO", version = "1.9"))
# Pkg.add(PackageSpec(name = "ForwardDiff", version = "0.10"))
# Pkg.add(PackageSpec(name = "IfElse", version = "0.1"))
# Pkg.add(PackageSpec(name = "JLD2", version = "0.4"))
# Pkg.add(PackageSpec(name = "Rotations", version = "1.0"))
# Pkg.add(PackageSpec(name = "StaticArrays", version = "1.2"))
# Pkg.add(PackageSpec(name = "Symbolics", version = "0.1.29"))
# Pkg.add(PackageSpec(name = "LinearAlgebra"))

# # Load minimal set of packages and scripts
# loaderdir = joinpath(@__DIR__, "..", "src", "loader.jl")
# include(loaderdir)

# # Generate dynamics expressions
# dynamicsdir = joinpath(@__DIR__, "..", "src", "dynamics")
# include(joinpath(dynamicsdir, "generate_dynamics.jl"))

# # Generate simulation expressions
# simulationdir = joinpath(@__DIR__, "..", "src", "simulation")
# include(joinpath(simulationdir, "generate_simulation2.jl"))
