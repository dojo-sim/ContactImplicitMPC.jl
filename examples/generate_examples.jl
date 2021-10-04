################################################################################
# Dynamics and simulation generation
################################################################################
# Generate dynamics
dynamicsdir = joinpath(@__DIR__, "..", "src", "dynamics")
include(joinpath(dynamicsdir, "generate_dynamics.jl"))

# Generate simulation
simulationdir = joinpath(@__DIR__, "..", "src", "simulation")
include(joinpath(simulationdir, "generate_simulation.jl"))
