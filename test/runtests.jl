using Test
using ContactControl
# using Colors
# using FFMPEG
using ForwardDiff
using JLD2
using MeshCat
using ModelingToolkit
using Parameters
using Plots
using Rotations
using StaticArrays
using LinearAlgebra
using Logging
using Random
using SparseArrays

# Dynamics
@testset "Dynamics Tests" begin
    include("dynamics/fast_model_methods.jl")
end

# # Simulator
# include("simulator/my_test_file.jl")
#
#
# # Controller
# include("controller/my_test_file.jl")
