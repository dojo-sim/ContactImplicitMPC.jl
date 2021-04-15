include("model.jl")
include("code_gen.jl")
include("fast_methods.jl")

################################################################################
# Particle
################################################################################
dir = joinpath(@__DIR__, "particle")
model = deepcopy(particle)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Particle 2D
################################################################################
dir = joinpath(@__DIR__, "particle_2D")
model = deepcopy(particle_2D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Hopper (2D)
################################################################################
dir = joinpath(@__DIR__, "hopper_2D")
model = deepcopy(hopper_2D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions_analytical(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Hopper (3D)
################################################################################
dir = joinpath(@__DIR__, "hopper_3D")
model = deepcopy(hopper_3D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions_analytical(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Quadruped
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Biped
################################################################################
dir = joinpath(@__DIR__, "biped")
model = deepcopy(biped)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Biped (5-link)
################################################################################
dir = joinpath(@__DIR__, "biped5")
model = deepcopy(biped5)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)
