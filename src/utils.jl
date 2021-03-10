# Run tests locally
using Pkg
Pkg.activate(@__DIR__)
# Pkg.activate(joinpath(@__DIR__, "test")

# Pkg.add("BenchmarkTools")
# Pkg.add("Colors")
# Pkg.add("FFMPEG")
# Pkg.add("ForwardDiff")
# Pkg.add("JLD2")
# Pkg.add("MeshCat")
# Pkg.add("ModelingToolkit")
# Pkg.add("Parameters")
# Pkg.add("Plots")
# Pkg.add("Rotations")
# Pkg.add("StaticArrays")
# Pkg.add("LinearAlgebra")
# Pkg.add("Logging")
# Pkg.add("Random")
# Pkg.add("SparseArrays")
# Pkg.add("Test")
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#
