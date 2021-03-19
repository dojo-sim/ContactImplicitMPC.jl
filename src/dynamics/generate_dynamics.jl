include("model.jl")
include("code_gen.jl")
include("fast_methods.jl")

################################################################################
# Particle (flat)
################################################################################
dir = joinpath(@__DIR__, "particle")
model = deepcopy(particle)

path_base = joinpath(dir, "flat/base.jld2")
path_dyn = joinpath(dir, "flat/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_approx = joinpath(dir, "flat/approx.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)

################################################################################
# Particle (quadratic)
################################################################################
dir = joinpath(@__DIR__, "particle")
model = deepcopy(particle_quadratic)

path_base = joinpath(dir, "quadratic/base.jld2")
path_dyn = joinpath(dir, "quadratic/dynamics.jld2")
path_res = joinpath(dir, "quadratic/residual.jld2")
path_jac = joinpath(dir, "quadratic/sparse_jacobians.jld2")
path_approx = joinpath(dir, "quadratic/approx.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)

################################################################################
# Particle 2D (flat)
################################################################################
dir = joinpath(@__DIR__, "particle_2D")
model = deepcopy(particle_2D)

path_base = joinpath(dir, "flat/base.jld2")
path_dyn = joinpath(dir, "flat/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_approx = joinpath(dir, "flat/approx.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)

################################################################################
# Particle 2D (slope)
################################################################################
dir = joinpath(@__DIR__, "particle_2D")
model = deepcopy(particle_2D_slope)
path_base = joinpath(dir, "slope/base.jld2")
path_dyn = joinpath(dir, "slope/dynamics.jld2")
path_res = joinpath(dir, "slope/residual.jld2")
path_jac = joinpath(dir, "slope/sparse_jacobians.jld2")
path_approx = joinpath(dir, "slope/approx.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)

################################################################################
# Hopper (2D)
################################################################################
dir = joinpath(@__DIR__, "hopper_2D")
model = deepcopy(hopper_2D)

path_base = joinpath(dir, "flat/base.jld2")
path_dyn = joinpath(dir, "flat/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_approx = joinpath(dir, "flat/approx.jld2")

expr_base = generate_base_expressions_analytical(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)

################################################################################
# Hopper (3D)
################################################################################
dir = joinpath(@__DIR__, "hopper_3D")
model = deepcopy(hopper_3D)

path_base = joinpath(dir, "flat/base.jld2")
path_dyn = joinpath(dir, "flat/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_approx = joinpath(dir, "flat/approx.jld2")

expr_base = generate_base_expressions_analytical(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)

################################################################################
# Quadruped
################################################################################
dir = joinpath(@__DIR__, "quadruped")
model = deepcopy(quadruped)
path_base = joinpath(dir, "flat/base.jld2")
path_dyn = joinpath(dir, "flat/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_approx = joinpath(dir, "flat/approx.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

expr_approx = generate_approx_expressions(model)
save_expressions(expr_approx, path_approx, overwrite=true)
instantiate_approx!(model, path_approx)
