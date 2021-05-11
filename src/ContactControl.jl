
module ContactControl

using BenchmarkTools
using Colors
using FFMPEG
using ForwardDiff
using GeometryBasics
using JLD2
using QDLDL
using MeshCat
using MeshCatMechanisms
using Meshing
using Symbolics
using IfElse
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
using FileIO
using MeshIO

# Utilities
include("utils.jl")

# Solver
include("solver/interior_point.jl")
include("solver/lu.jl")
include("solver/ldl.jl")
include("solver/qr.jl")
include("solver/schur.jl")

# Environment
include("simulator/environment.jl")

# Dynamics
include("dynamics/model.jl")
include("dynamics/code_gen.jl")
include("dynamics/fast_methods.jl")

# Simulator
include("simulator/trajectory.jl")

export ContactDynamicsModel, Dimensions, BaseMethods, DynamicsMethods,
    ResidualMethods, Environment
export environment_2D, environment_3D, environment_2D_flat,
    environment_3D_flat, get_model

# Models
include("dynamics/particle_2D/model.jl")
include("dynamics/particle/model.jl")
include("dynamics/hopper_2D/model.jl")
include("dynamics/hopper_3D/model.jl")
include("dynamics/quadruped/model.jl")
include("dynamics/biped/model.jl")
include("dynamics/biped5/model.jl")
include("dynamics/flamingo/model.jl")
include("dynamics/pushbot/model.jl")

# Simulator
include("simulator/policy.jl")
include("simulator/disturbances.jl")
include("simulator/simulator.jl")

# Controller
include("controller/linearized_step.jl")
include("controller/implicit_dynamics.jl")
include("controller/objective.jl")
include("controller/linearized_solver.jl")
include("controller/newton.jl")
include("controller/mpc_utils.jl")
include("controller/policy.jl")

# Visuals
include("dynamics/visuals.jl")
include("dynamics/visual_utils.jl")

export SparseStructure, LinearizedStep, get_bilinear_indices,
    bil_addition!, r_linearized!, rz_linearized!

export ImplicitTraj, linearization!, implicit_dynamics!
export TrackingObjective

end # module
