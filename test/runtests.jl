using Test
using ContactControl
using ForwardDiff
using JLD2
using ModelingToolkit
using StaticArrays
using LinearAlgebra
using Random
using SparseArrays

# Dynamics
include("dynamics/particle.jl")
# include("dynamics/quadruped.jl")

# Simulator
include("simulator/random_qp.jl")
include("simulator/particle.jl")

# Controller
include("controller/bilinear.jl")
include("controller/implicit_dynamics.jl")
include("controller/newton.jl")

# Gram Schmidt
include("cgs/cgs.jl")
