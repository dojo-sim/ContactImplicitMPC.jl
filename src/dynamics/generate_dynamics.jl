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

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
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

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
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

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Quadruped Payload
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped_payload)
path_base = joinpath(dir, "dynamics_payload/base.jld2")
path_dyn = joinpath(dir, "dynamics_payload/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Flamingo
################################################################################
dir = joinpath(@__DIR__, "flamingo")
model = deepcopy(flamingo)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# PushBot
################################################################################
dir = joinpath(@__DIR__, "pushbot")
model = deepcopy(pushbot)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Walled Cartpole
################################################################################
dir = joinpath(@__DIR__, "walledcartpole")
model = deepcopy(walledcartpole)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=false, C_analytical=false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

###############################################################################
# Rigid body
###############################################################################
dir = joinpath(@__DIR__, "rigidbody")
model = deepcopy(rigidbody)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model,
	M_analytical = true,
	C_analytical = true,
	mapping = G_func,
	nv = model.nq - 1)

save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model, nv = model.nq - 1)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Centroidal Quadruped (3D)
################################################################################
dir = joinpath(@__DIR__, "centroidal_quadruped")
model = deepcopy(centroidal_quadruped)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Centroidal Quadruped Wall (3D)
################################################################################
dir = joinpath(@__DIR__, "centroidal_quadruped_wall")
model = deepcopy(centroidal_quadruped_wall)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Centroidal Quadruped (3D) undamped
################################################################################
dir = joinpath(@__DIR__, "centroidal_quadruped")
model = deepcopy(centroidal_quadruped_undamped)
path_base = joinpath(dir, "dynamics_undamped/base.jld2")
path_dyn = joinpath(dir, "dynamics_undamped/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

################################################################################
# Point Foot Quadruped
################################################################################
dir = joinpath(@__DIR__, "point_foot_quadruped")
model = deepcopy(point_foot_quadruped)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

# ################################################################################
# # Point Foot Quadruped torsional
# ################################################################################
# dir = joinpath(@__DIR__, "point_foot_quadruped")
# model = deepcopy(point_foot_quadruped)
# path_base = joinpath(dir, "dynamics/base.jld2")
# path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
#
# expr_base = generate_base_expressions(model, M_analytical=true, C_analytical=true)
# save_expressions(expr_base, path_base, overwrite=true)
# instantiate_base!(model, path_base)
#
# expr_dyn = generate_dynamics_expressions(model)
# save_expressions(expr_dyn, path_dyn, overwrite=true)
# instantiate_dynamics!(model, path_dyn)
