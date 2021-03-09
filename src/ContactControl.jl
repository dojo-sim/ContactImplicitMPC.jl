
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

# Models
include("dynamics/particle/model.jl")
# include("dynamics/quadruped/model.jl")

export ContactDynamicsModel, Dimensions, BaseMethods, DynamicsMethods, ResidualMethods

# Simulator
include("simulator/interior_point.jl")
include("simulator/simulator.jl")

# COntroller
include("controller/bilinear.jl")

end # module
