
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

include(joinpath(pwd(), "src/dynamics/model.jl"))
include(joinpath(pwd(), "src/dynamics/code_gen.jl"))

end # module
