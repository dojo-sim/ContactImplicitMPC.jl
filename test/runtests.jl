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
include("dynamics/fast_model_methods.jl")

# Simulator
include("simulator/random_qp.jl")
include("simulator/particle.jl")

# Controller
include("controller/bilinear.jl")
include("controller/implicit_dynamics.jl")
# include("controller/newton.jl")
