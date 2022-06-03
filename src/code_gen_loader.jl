include(joinpath(@__DIR__, ".."))

using FileIO
using ForwardDiff
using IfElse
using JLD2
using LinearAlgebra
using Rotations
using RoboDojo
import RoboDojo: LinearSolver, LUSolver, Model, ResidualMethods, Space, Disturbances, IndicesZ, InteriorPoint, EmptySolver, Policy, Trajectory, GradientTrajectory, InteriorPointOptions, IndicesOptimization, interior_point, interior_point_solve!, bilinear_violation, residual_violation, general_correction_term!, r!, rz!, rθ!, linear_solve!, lu_solver, empty_policy, empty_disturbances, friction_coefficients, SimulatorStatistics, SimulatorOptions, indices_θ, num_data, initialize_z!, initialize_θ!, indices_z, indices_θ, simulate!, policy, process!, Simulator

using StaticArrays
using Symbolics

# Utilities
include("utils.jl")

# Solver

# Environment
include("simulator/environment.jl")

# Dynamics
include("dynamics/model.jl")

# Simulator
include("simulation/index.jl")

# Simulator
include("simulation/contact_methods.jl")
include("simulation/simulation.jl")

include("dynamics/code_gen_dynamics.jl")
include("dynamics/fast_methods_dynamics.jl")

# # Models
include("dynamics/particle_2D/model.jl")
include("dynamics/particle/model.jl")
include("dynamics/hopper_2D/model.jl")
include("dynamics/hopper_3D/model.jl")
include("dynamics/quadruped/model.jl")
include("dynamics/flamingo/model.jl")
include("dynamics/pushbot/model.jl")
include("dynamics/walledcartpole/model.jl")
include("dynamics/centroidal_quadruped/model.jl")
include("dynamics/centroidal_quadruped_wall/model.jl")
include("dynamics/centroidal_quadruped_box/model.jl")

# Simulation
include("simulation/environments/flat.jl")
include("simulation/environments/piecewise.jl")
include("simulation/environments/quadratic.jl")
include("simulation/environments/slope.jl")
include("simulation/environments/sinusoidal.jl")
include("simulation/environments/stairs.jl")

include("simulation/residual_approx.jl")
include("simulation/code_gen_simulation.jl")
