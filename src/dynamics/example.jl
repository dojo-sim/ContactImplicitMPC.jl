include("model_struct.jl")
include("quadruped_model.jl")
include("code_gen.jl")

T = Float64

model = deepcopy(quadruped)
# expr = generate_dynamics_expressions(model)
save_dynamics_expressions(expr, joinpath(@__DIR__, ".expr/quadruped_expr.jld2"), overwrite=true)
load_dynamics_expressions(joinpath(@__DIR__, ".expr/quadruped_expr.jld2"))
instantiate_dynamics!(model, joinpath(@__DIR__, ".expr/quadruped_expr.jld2"))


nq = model.dim.q
nu = model.dim.u
nγ = model.dim.γ
nb = model.dim.b
nz = model.dim.z

q0s = rand(SizedVector{nq,T})
q1s = rand(SizedVector{nq,T})
u1s = rand(SizedVector{nu,T})
γ1s = rand(SizedVector{nγ,T})
b1s = rand(SizedVector{nb,T})
q2s = rand(SizedVector{nq,T})
q̇0s = rand(SizedVector{nq,T})

M1s   = rand(SizedMatrix{nq,nq,T})
∇zs   = rand(SizedMatrix{nq,nz,T})
∇q0s  = rand(SizedMatrix{nq,nq,T})
∇q1s  = rand(SizedMatrix{nq,nq,T})
∇u1s  = rand(SizedMatrix{nq,nu,T})
∇γ1s  = rand(SizedMatrix{nq,nγ,T})
∇b1s  = rand(SizedMatrix{nq,nb,T})
∇q2s  = rand(SizedMatrix{nq,nq,T})



∇z_dynamics_fast!(model, ∇zs, q0s, q1s, u1s, γ1s, b1s, q2s)
∇q_1_dynamics_fast!(model, ∇q0s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇q_dynamics_fast!(model,   ∇q1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇u_dynamics_fast!(model,   ∇u1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇γ_dynamics_fast!(model,   ∇γ1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇b_dynamics_fast!(model,   ∇b1s, q0s, q1s, u1s, γ1s, b1s, q2s)
∇q1_dynamics_fast!(model,  ∇q2s, q0s, q1s, u1s, γ1s, b1s, q2s)

∇zs
∇q0s
∇q1s
∇u1s
∇γ1s
∇b1s
∇q2s







quadruped
instantiate_dynamics!(quadruped, joinpath(@__DIR__, "quadruped_expr.jld2"))
lagrangian_fast(quadruped, q1s, q̇1s)
M_fast(quadruped, q1s)
B_fast(quadruped, q1s)
N_fast(quadruped, q1s)
P_fast(quadruped, q1s)
C_fast(quadruped, q1s, q̇1s)
dynamics_fast(quadruped, q1s, q2s, u2s, γ2s, b2s, q3s)
∇z_dynamics_fast!(quadruped,   ∇zs, q1s, q2s, u2s, γ2s, b2s, q3s)
∇q_1_dynamics_fast!(quadruped, ∇q_1s, q1s, q2s, u2s, γ2s, b2s, q3s)
∇q_dynamics_fast!(quadruped,   ∇qs, q1s, q2s, u2s, γ2s, b2s, q3s)
∇u_dynamics_fast!(quadruped,   ∇us, q1s, q2s, u2s, γ2s, b2s, q3s)
∇γ_dynamics_fast!(quadruped,   ∇γs, q1s, q2s, u2s, γ2s, b2s, q3s)
∇b_dynamics_fast!(quadruped,   ∇bs, q1s, q2s, u2s, γ2s, b2s, q3s)
∇q1_dynamics_fast!(quadruped,  ∇qs, q1s, q2s, u2s, γ2s, b2s, q3s)
