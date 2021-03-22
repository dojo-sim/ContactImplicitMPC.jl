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

# Solver
include("solver/gs.jl")
include("solver/qdldl.jl")
include("solver/random_qp.jl")

# Dynamics
include("dynamics/particle.jl")
include("dynamics/quadruped.jl")

# Simulator
include("simulator/particle.jl")
include("simulator/hopper_2D.jl")
include("simulator/hopper_3D.jl")


# Controller
include("controller/cost_function.jl")
include("controller/bilinear.jl")
include("controller/implicit_dynamics.jl")
# include("controller/newton.jl")
