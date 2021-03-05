include("model.jl")
include("quadruped/model.jl")
include("code_gen.jl")

L_fast, M_fast, B_fast, N_fast, P_fast, C_fast, d, dy!, dq0!, dq1!, du1!, dγ1!, db1!, dq2!, r!, rz!, rθ!, rz_sp, rθ_sp = generate_fast_expressions(model, "quadruped", generate = false);

# nq = model.dim.q
# nu = model.dim.u
# nc = model.dim.c
# nb = model.dim.b
#
# q0s = rand(SizedVector{nq,T})
# q1s = rand(SizedVector{nq,T})
# u1s = rand(SizedVector{nu,T})
# γ1s = rand(SizedVector{nγ,T})
# b1s = rand(SizedVector{nb,T})
# q2s = rand(SizedVector{nq,T})
# q̇0s = rand(SizedVector{nq,T})
#
# M1s   = rand(SizedMatrix{nq,nq,T})
# ∇ys   = rand(SizedMatrix{nq,ny,T})
# ∇q0s  = rand(SizedMatrix{nq,nq,T})
# ∇q1s  = rand(SizedMatrix{nq,nq,T})
# ∇u1s  = rand(SizedMatrix{nq,nu,T})
# ∇γ1s  = rand(SizedMatrix{nq,nγ,T})
# ∇b1s  = rand(SizedMatrix{nq,nb,T})
# ∇q2s  = rand(SizedMatrix{nq,nq,T})
#
#
# M_fast(model, q0s)
# B_fast(model, q0s)
# N_fast(model, q0s)
# P_fast(model, q0s)
# C_fast(model, q0s, q̇0s)
# dynamics_fast(model, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇y_dynamics_fast!(model, ∇ys, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇q0_dynamics_fast!(model, ∇q0s, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇q1_dynamics_fast!(model,   ∇q1s, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇u1_dynamics_fast!(model,   ∇u1s, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇γ1_dynamics_fast!(model,   ∇γ1s, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇b1_dynamics_fast!(model,   ∇b1s, q0s, q1s, u1s, γ1s, b1s, q2s)
# ∇q2_dynamics_fast!(model,  ∇q2s, q0s, q1s, u1s, γ1s, b1s, q2s)
#
# nz = model.dim.z
# nθ = model.dim.θ
# zs = rand(SizedVector{nz})
# θs = rand(SizedVector{nθ})
# κs = 1e-3
# rs = rand(SizedVector{nz})
# ∇zs = rand(nz,nz)
# ∇θs = rand(nz,nθ)
#
#
# r_fast!(model, rs, zs, θs, κs)
# rz_fast!(model, ∇zs, zs, θs, κs)
# rθ_fast!(model, ∇θs, zs, θs, κs)
#
#
# function line_search_eval3!(model::ContactDynamicsModel, out::Vector{SizedArray{Tuple{nz},T,1,1}},
#     z::SizedVector{nz,T}, δz::SizedVector{nz,T}, θ::SizedVector{nθ,T}, κ::T, N::Int) where {nz,nθ,T}
#     for k = 1:N
#         r_fast!(model, out[k], z+1/2^k*δz, θ, κ)
#     end
#     return nothing
# end
#
# function multi_line_search_eval3!(model::ContactDynamicsModel, out::Vector{SizedArray{Tuple{nz},T,1,1}},
#     z::SizedVector{nz,T}, δz::SizedVector{nz,T}, θ::SizedVector{nθ,T}, κ::T, N::Int) where {nz,nθ,T}
#     Threads.@threads for k = 1:N
#         r_fast!(model, out[k], z+1/2^k*δz, θ, κ)
#     end
#     return nothing
# end
#
#
# model = deepcopy(quadruped)
# instantiate_base!(model, joinpath(@__DIR__, ".expr/quadruped_base.jld2"))
# instantiate_dynamics!(model, joinpath(@__DIR__, ".expr/quadruped_dynamics.jld2"))
# instantiate_residual!(model, joinpath(@__DIR__, ".expr/quadruped_residual.jld2"))
#
# N = 12
# outs0 = [rand(SizedVector{nz,T}) for k=1:N]
# outs1 = [rand(SizedVector{nz,T}) for k=1:N]
# zs = rand(SizedVector{nz,T})
# δzs = rand(SizedVector{nz,T})
# θs = rand(SizedVector{nθ,T})
# κs = 1e-3
# @time line_search_eval3!(model, outs0, zs, δzs, θs, κs, N)
# @time multi_line_search_eval3!(model, outs1, zs, δzs, θs, κs, N)
#
# @btime line_search_eval3!(model, outs0, zs, δzs, θs, κs, N)
# @btime multi_line_search_eval3!(model, outs1, zs, δzs, θs, κs, N)
#
# outs0 == outs1
#
# rz_fast!(model, ∇zs, zs, θs, κs)
# @btime ∇zs\rs
