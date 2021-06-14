
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

using Ipopt, MathOptInterface


# Utilities
include("utils.jl")

# Solver
include("solver/cones.jl")
include("solver/interior_point.jl")
include("solver/lu.jl")
include("solver/gn.jl")
include("solver/ldl.jl")
include("solver/qr.jl")
include("solver/schur.jl")

# Environment
include("simulator/environment.jl")

# Dynamics
include("dynamics/model.jl")

# Simulator
include("simulation/contact_methods.jl")
include("simulation/simulation.jl")
include("simulator/trajectory.jl")

include("dynamics/code_gen_dynamics.jl")
include("dynamics/fast_methods_dynamics.jl")

export ContactModel, Dimensions, BaseMethods, DynamicsMethods,
    ResidualMethods, Environment
export environment_2D, environment_3D, environment_2D_flat,
    environment_3D_flat, get_model

# Models
include("dynamics/quaternions.jl")
include("dynamics/mrp.jl")
include("dynamics/euler.jl")

include("dynamics/particle_2D/model.jl")
include("dynamics/particle/model.jl")
include("dynamics/hopper_2D/model.jl")
include("dynamics/hopper_3D/model.jl")
include("dynamics/hopper_3D_quaternion/model.jl")
include("dynamics/quadruped/model.jl")
include("dynamics/quadruped_simple/model.jl")
include("dynamics/biped/model.jl")
include("dynamics/flamingo/model.jl")
include("dynamics/pushbot/model.jl")
include("dynamics/planarpush/model.jl")
include("dynamics/rigidbody/model.jl")
include("dynamics/box/model.jl")

# Simulator
include("simulator/policy.jl")
include("simulator/disturbances.jl")
include("simulator/simulator.jl")

# Simulation
include("simulation/environments/flat.jl")
include("simulation/environments/piecewise.jl")
include("simulation/environments/quadratic.jl")
include("simulation/environments/slope.jl")
include("simulation/environments/sinusoidal.jl")
include("simulation/environments/stairs.jl")

include("simulation/residual_approx.jl")
include("simulation/code_gen_simulation.jl")

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
