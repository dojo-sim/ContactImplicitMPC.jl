
module ContactControl

greet() = print("ContactControl")

using BenchmarkTools
using Colors
using FFMPEG
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
using Test

# Dynamics
include("dynamics/model.jl")
include("dynamics/code_gen.jl")
include("dynamics/fast_methods.jl")

export ContactDynamicsModel, Dimensions, BaseMethods, DynamicsMethods, ResidualMethods

# Simulator
include("simulator/interior_point.jl")
include("simulator/simulator.jl")

# Controller
include("controller/bilinear.jl")

export SparseStructure, LinStep, get_bilinear_indices, bil_addition!, r_approx!, rz_approx!

# Models
include("dynamics/particle/model.jl")
# include("dynamics/quadruped/model.jl")

end # module
