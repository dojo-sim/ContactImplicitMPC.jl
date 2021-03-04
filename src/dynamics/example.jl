include("model_struct.jl")
include("quadruped_model.jl")
include("code_gen.jl")
include("model_methods.jl")
include("fast_model_methods.jl")

T = Float64

model = deepcopy(quadruped)
# expr_bas = generate_base_expressions(model)
# save_expressions(expr_bas, joinpath(@__DIR__, ".expr/quadruped_base.jld2"), overwrite=true)
instantiate_base!(model, joinpath(@__DIR__, ".expr/quadruped_base.jld2"))

# expr_dyn = generate_dynamics_expressions(model)
# save_expressions(expr_dyn, joinpath(@__DIR__, ".expr/quadruped_dynamics.jld2"), overwrite=true)
instantiate_dynamics!(model, joinpath(@__DIR__, ".expr/quadruped_dynamics.jld2"))

# expr_res = generate_residual_expressions(model)
# save_expressions(expr_res, joinpath(@__DIR__, ".expr/quadruped_residual.jld2"), overwrite=true)
instantiate_residual!(model, joinpath(@__DIR__, ".expr/quadruped_residual.jld2"))


nq = model.dim.q
nu = model.dim.u
nγ = model.dim.γ
nb = model.dim.b
ny = model.dim.y
nz = model.dim.z
nθ = model.dim.θ

q0s = rand(SizedVector{nq,T})
q1s = rand(SizedVector{nq,T})
u1s = rand(SizedVector{nu,T})
γ1s = rand(SizedVector{nγ,T})
b1s = rand(SizedVector{nb,T})
q2s = rand(SizedVector{nq,T})
q̇0s = rand(SizedVector{nq,T})

M1s   = rand(SizedMatrix{nq,nq,T})
∇ys   = rand(SizedMatrix{nq,ny,T})
∇q0s  = rand(SizedMatrix{nq,nq,T})
∇q1s  = rand(SizedMatrix{nq,nq,T})
∇u1s  = rand(SizedMatrix{nq,nu,T})
∇γ1s  = rand(SizedMatrix{nq,nγ,T})
∇b1s  = rand(SizedMatrix{nq,nb,T})
∇q2s  = rand(SizedMatrix{nq,nq,T})


M_fast(model, q0s)
B_fast(model, q0s)
N_fast(model, q0s)
P_fast(model, q0s)
C_fast(model, q0s, q̇0s)
dynamics_fast(model, q0s, q1s, u1s, γ1s, b1s, q2s)
∇y_dynamics_fast!(model, ∇ys, q0s, q1s, u1s, γ1s, b1s, q2s)
∇q0_dynamics_fast!(model, ∇q0s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇q1_dynamics_fast!(model,   ∇q1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇u1_dynamics_fast!(model,   ∇u1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇γ1_dynamics_fast!(model,   ∇γ1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇b1_dynamics_fast!(model,   ∇b1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇q2_dynamics_fast!(model,  ∇q2s, q0s, q1s, u1s, γ1s, b1s, q2s)




model = deepcopy(quadruped)
instantiate_base!(model, joinpath(@__DIR__, ".expr/quadruped_base.jld2"))
instantiate_dynamics!(model, joinpath(@__DIR__, ".expr/quadruped_dynamics.jld2"))
# expr = generate_residual_expressions(model)
# save_expressions(expr, joinpath(@__DIR__, ".expr/quadruped_residual.jld2"), overwrite=true)
instantiate_residual!(model, joinpath(@__DIR__, ".expr/quadruped_residual.jld2"))


zs = rand(SizedVector{nz})
θs = rand(SizedVector{nθ})
ρs = 1e-3
rs = rand(SizedVector{nz})
∇zs = rand(nz,nz)
∇θs = rand(nz,nθ)
@btime r_fast!(model, rs, zs, θs, ρs)
@btime rz_fast!(model, ∇zs, zs, θs, ρs)
@btime rθ_fast!(model, ∇θs, zs, θs, ρs)
