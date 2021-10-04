include(joinpath(@__DIR__, ".."))

using FileIO
using ForwardDiff
using IfElse
using JLD2
using LinearAlgebra
using Rotations
using StaticArrays
using Symbolics

# Utilities
include("utils.jl")

# Solver
include("solver/cones.jl")

# Environment
include("simulator/environment.jl")

# Dynamics
include("dynamics/model.jl")

# Simulator
include("simulation/index.jl")

# Simulator
include("simulation/contact_methods.jl")
include("simulation/simulation.jl")
include("simulator/trajectory.jl")

include("dynamics/code_gen_dynamics.jl")
include("dynamics/fast_methods_dynamics.jl")

# # Models
include("dynamics/quaternions.jl")

include("dynamics/particle_2D/model.jl")
include("dynamics/particle/model.jl")
include("dynamics/hopper_2D/model.jl")
include("dynamics/hopper_3D/model.jl")
include("dynamics/quadruped/model.jl")
include("dynamics/flamingo/model.jl")
include("dynamics/pushbot/model.jl")
include("dynamics/rigidbody/model.jl")

# Simulation
include("simulation/environments/flat.jl")
include("simulation/environments/piecewise.jl")
include("simulation/environments/quadratic.jl")
include("simulation/environments/slope.jl")
include("simulation/environments/sinusoidal.jl")
include("simulation/environments/stairs.jl")

include("simulation/residual_approx.jl")
include("simulation/code_gen_simulation.jl")
