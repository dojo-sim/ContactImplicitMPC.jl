using Test
using ContactControl
# using Colors
# # using FFMPEG
using ForwardDiff
using JLD2
# using MeshCat
using ModelingToolkit
# using Parameters
# using Plots
# using Rotations
using StaticArrays
using LinearAlgebra
# using Logging
using Random
using SparseArrays

# Dynamics
@testset "Dynamics Tests" begin
    include("dynamics/fast_model_methods.jl")
end

# Simulator
@testset "Simulator Tests" begin
    include("simulator/random_qp.jl")
    include("simulator/particle.jl")
end

# # Controller
# include("controller/my_test_file.jl")
