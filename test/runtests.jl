using Test
using ForwardDiff
using JLD2
using Symbolics
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using ContactImplicitMPC

# Solver
include("solver/qr.jl")
include("solver/lu.jl")
include("solver/schur.jl")

# Dynamics
include("dynamics/model.jl")
include("dynamics/particle.jl")
include("dynamics/quadruped.jl")

# Simulator
include("simulator/rotations.jl")
include("simulator/environment.jl")
include("simulator/particle.jl")
include("simulator/hopper_3D.jl")
include("simulator/quadruped.jl")
include("simulator/open_loop.jl")

# Controller
include("controller/objective.jl")
include("controller/linearized_step.jl")
include("controller/implicit_dynamics.jl")
include("controller/linearized_solver.jl")
include("controller/newton.jl")
include("controller/newton_structure_solver.jl")
include("controller/trajectory.jl")

# MPC examples
include("controller/mpc_quadruped.jl")
include("controller/mpc_flamingo.jl") #TODO fix newton structure solver test
