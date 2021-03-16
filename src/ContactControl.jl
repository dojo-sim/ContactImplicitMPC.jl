
module ContactControl

greet() = print("ContactControl")

using BenchmarkTools
using Colors
using FFMPEG
using ForwardDiff
using GeometryBasics
using JLD2
using QDLDL
using MeshCat
using ModelingToolkit
using Parameters
using Plots
using Rotations
using CoordinateTransformations
using StaticArrays
using LinearAlgebra
using Logging
using Random
using SparseArrays
using Test

# Solver
include("solver/interior_point.jl")

# Dynamics
include("dynamics/environment.jl")
include("dynamics/model.jl")
include("dynamics/code_gen.jl")
include("dynamics/fast_methods.jl")
include("dynamics/visuals.jl")

export ContactDynamicsModel, Dimensions, BaseMethods, DynamicsMethods, ResidualMethods, Environment
export environment_2D, environment_3D, environment_2D_flat, environment_3D_flat, get_model

# Models
include("dynamics/particle_2D/model.jl")
include("dynamics/particle/model.jl")
include("dynamics/quadruped/model.jl")

# Simulator
include("simulator/trajectory.jl")
include("simulator/simulator.jl")

# Controller
include("controller/bilinear.jl")
include("controller/implicit_dynamics.jl")
include("controller/cost_function.jl")
include("controller/qdldl.jl")
# include("controller/newton.jl")

export SparseStructure, LinStep, get_bilinear_indices, bil_addition!, r_approx!, rz_approx!
export ImplicitTraj, linearization!, implicit_dynamics!
export CostFunction

# CGS
# include("gram_schmidt/gram_schmidt.jl")

# export triangularize, triu_perm, CGSData13, CGSData13!, cgs!, qr_solve!

end # module
