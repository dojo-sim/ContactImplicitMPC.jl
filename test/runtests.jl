using Test
using ContactControl
using ForwardDiff
using JLD2
using QDLDL
using Symbolics
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using Quaternions

Solver
include("solver/gs.jl")
include("solver/qdldl.jl")
include("solver/random_qp.jl")
include("solver/lu.jl")
include("solver/gs.jl")
include("solver/schur.jl")
include("solver/soc.jl")
include("solver/mehrotra.jl")

# Dynamics
# include("dynamics/lagrangian.jl") #NEED FIX
include("dynamics/model.jl")
include("dynamics/particle.jl")
include("dynamics/quadruped.jl")
include("dynamics/quaternion.jl")

# Simulator
include("simulator/rotations.jl")
include("simulator/environment.jl")
include("simulator/trajectory.jl")
# include("simulator/particle.jl") #NEED FIX
# include("simulator/hopper_2D.jl") #TODO: set tests to raibert model #NEED FIX
# include("simulator/hopper_3D.jl") #NEED FIX
include("simulator/quadruped.jl")
include("simulator/open_loop.jl")
# include("simulator/biped.jl") #TODO: improve this test #NEED FIX
# include("simulator/flamingo.jl") #TODO: add this test #NEED FIX

# Controller
include("controller/objective.jl")
include("controller/linearized_step.jl")
include("controller/implicit_dynamics.jl")
# include("controller/linearized_solver.jl") #NEED FIX
include("controller/newton.jl")

# MPC examples
include("controller/mpc_quadruped.jl")
include("controller/mpc_flamingo.jl")
# const ContactControl = Main
