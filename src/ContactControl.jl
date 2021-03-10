
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
include("dynamics/environment.jl")
include("dynamics/model.jl")
include("dynamics/code_gen.jl")
include("dynamics/fast_methods.jl")
include("dynamics/visuals.jl")

# Models
include("dynamics/particle/model.jl")
include("dynamics/quadruped/model.jl")

export ContactDynamicsModel, Dimensions, BaseMethods, DynamicsMethods, ResidualMethods, Environment
export environment_2D, environment_3D, environment_2D_flat, environment_3D_flat

# Simulator
include("simulator/interior_point.jl")
include("simulator/simulator.jl")

end # module
