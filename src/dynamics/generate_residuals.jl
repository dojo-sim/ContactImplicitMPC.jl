include("model.jl")
include("code_gen.jl")
include("fast_methods.jl")

################################################################################
# Particle (flat)
################################################################################
dir = joinpath(@__DIR__, "particle")
model = deepcopy(particle)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Particle (quadratic)
################################################################################
dir = joinpath(@__DIR__, "particle")
model = deepcopy(particle_quadratic)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "quadratic/residual.jld2")
path_jac = joinpath(dir, "quadratic/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "quadratic/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Particle (sinusoidal)
################################################################################
dir = joinpath(@__DIR__, "particle")
model = deepcopy(particle_sinusoidal)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Particle 2D (flat)
################################################################################
dir = joinpath(@__DIR__, "particle_2D")
model = deepcopy(particle_2D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Particle 2D (slope)
################################################################################
dir = joinpath(@__DIR__, "particle_2D")
model = deepcopy(particle_2D_slope)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "slope/residual.jld2")
path_jac = joinpath(dir, "slope/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "slope/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Hopper (2D)
################################################################################
dir = joinpath(@__DIR__, "hopper_2D")
model = deepcopy(hopper_2D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Hopper (2D) vertical hop
################################################################################
dir = joinpath(@__DIR__, "hopper_2D")
model = deepcopy(hopper_2D_vertical)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "vertical/residual.jld2")
path_jac = joinpath(dir, "vertical/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "vertical/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Hopper (2D) (sinusoidal)
################################################################################
dir = joinpath(@__DIR__, "hopper_2D")
model = deepcopy(hopper_2D_sinusoidal)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Hopper (2D) (stairs)
################################################################################
dir = joinpath(@__DIR__, "hopper_2D")
model = deepcopy(hopper_2D_stairs)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "stairs/residual.jld2")
path_jac = joinpath(dir, "stairs/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "stairs/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Hopper (3D)
################################################################################
dir = joinpath(@__DIR__, "hopper_3D")
model = deepcopy(hopper_3D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Hopper (3D) (sinusoidal)
################################################################################
dir = joinpath(@__DIR__, "hopper_3D")
model = deepcopy(hopper_3D_sinusoidal)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Quadruped
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Quadruped Payload
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped_payload)

path_base = joinpath(dir, "dynamics_payload/base.jld2")
path_dyn = joinpath(dir, "dynamics_payload/dynamics.jld2")
path_res = joinpath(dir, "payload/residual.jld2")
path_jac = joinpath(dir, "payload/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "payload/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Quadruped Sinusoidal
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped_sinusoidal)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Quadruped (Piecewise)
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped_piecewise)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "piecewise/residual.jld2")
path_jac = joinpath(dir, "piecewise/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "piecewise/linearized.jld2")

instantiate_base!(model_sim, path_base)

expr_dyn = generate_dynamics_expressions(model_sim, derivs = true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model_sim, path_dyn, derivs = true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model_sim, jacobians = :approx)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model_sim, path_res, jacobians = :approx)

################################################################################
# Quadruped (3D)
################################################################################
dir = joinpath(@__DIR__, "quadruped_3D")
model = deepcopy(quadruped_3D)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Quadruped (Linear)
################################################################################
dir = joinpath(@__DIR__, "quadrupedlinear")
model = deepcopy(quadrupedlinear)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)
model.spa.rz_sp = rz_sp
model.spa.rθ_sp = rθ_sp

################################################################################
# Biped
################################################################################
dir = joinpath(@__DIR__, "biped")
model = deepcopy(biped)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Biped (sinusoidal)
################################################################################
dir = joinpath(@__DIR__, "biped")
model = deepcopy(biped_sinusoidal)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Flamingo (flat)
################################################################################
dir = joinpath(@__DIR__, "flamingo")
model = deepcopy(flamingo)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Flamingo (sinusoidal)
################################################################################
dir = joinpath(@__DIR__, "flamingo")
model = deepcopy(flamingo_sinusoidal)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# Flamingo (smooth slope)
################################################################################
dir = joinpath(@__DIR__, "flamingo")
model = deepcopy(flamingo_smooth_slope)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "smooth_slope/residual.jld2")
path_jac = joinpath(dir, "smooth_slope/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "smooth_slope/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# PushBot
################################################################################
dir = joinpath(@__DIR__, "pushbot")
model = deepcopy(pushbot)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

################################################################################
# PlanarPush
################################################################################
dir = joinpath(@__DIR__, "planarpush")
model = deepcopy(planarpush)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)
