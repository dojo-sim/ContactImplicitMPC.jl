################################################################################
# Dynamics and simulation generation
################################################################################
using Pkg
pkg_dir = joinpath(dirname(pathof(ContactImplicitMPC)), "..")
Pkg.activate(pkg_dir)
include(joinpath(pkg_dir, "src/code_gen_loader.jl"))

# Generate dynamics
dynamicsdir = joinpath(dirname(pathof(ContactImplicitMPC)), "..", "src", "dynamics")
include(joinpath(dynamicsdir, "generate_dynamics.jl"))

# Generate simulation
simulationdir = joinpath(dirname(pathof(ContactImplicitMPC)), "..", "src", "simulation")
include(joinpath(simulationdir, "generate_simulation.jl"))