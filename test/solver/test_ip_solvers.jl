
# Test basic IP solver and Mehrotra base solver
include(joinpath(module_dir(), "src", "solver", "interior_point.jl"))
include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))
include("mehrotra.jl")

# Test Mehrotra expanded solver
include(joinpath(module_dir(), "src", "solver", "mehrotra_expanded.jl"))
include("mehrotra_expanded.jl")

# Benchmark Mehrotra base


# Benchmark Mehrotra expanded
