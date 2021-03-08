include("model.jl")
include("code_gen.jl")
include("fast_methods.jl")

################################################################################
# Quadruped
################################################################################
include("quadruped/model.jl")
model = deepcopy(quadruped)

dir = jointpath(@__DIR__, "quadruped")

path_base = joinpath(dir, "base.jld2")
path_dyn = joinpath(dir, "dynamics.jld2")
path_res = joinpath(dir, "residual.jld2")
path_jac = joinpath(dir, "sparse_jacobians.jld2")

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



################################################################################
# Particle
################################################################################
include("particle/model.jl")
model = deepcopy(particle)

dir = jointpath(@__DIR__, "particle")

path_base = joinpath(dir, "base.jld2")
path_dyn = joinpath(dir, "dynamics.jld2")
path_res = joinpath(dir, "residual.jld2")
path_jac = joinpath(dir, "sparse_jacobians.jld2")

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
