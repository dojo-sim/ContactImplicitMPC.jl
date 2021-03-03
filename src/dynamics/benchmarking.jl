include("quadruped_model.jl")

Pkg.add("StaticArrays")
Pkg.add("ForwardDiff")
Pkg.add("BenchmarkTools")
Pkg.add("JLD2")

using StaticArrays
using ForwardDiff
using ModelingToolkit
using BenchmarkTools
using LinearAlgebra
using JLD2


T = Float64

model = quadruped
nq = model.dim.q
nu = model.dim.u
nγ = model.dim.γ
nb = model.dim.b

q1 = rand(SVector{nq,T})
q2 = rand(SVector{nq,T})
u2 = rand(SVector{nu,T})
γ2 = rand(SVector{nγ,T})
b2 = rand(SVector{nb,T})
q3 = rand(SVector{nq,T})
q̇1 = rand(SVector{nq,T})

q1s = rand(SizedVector{nq,T})
q2s = rand(SizedVector{nq,T})
u2s = rand(SizedVector{nu,T})
γ2s = rand(SizedVector{nγ,T})
b2s = rand(SizedVector{nb,T})
q3s = rand(SizedVector{nq,T})
q̇1s = rand(SizedVector{nq,T})

M1s   = rand(SizedMatrix{nq,nq,T})
∇zs   = rand(SizedMatrix{nq,nq,T})
∇q_1s = rand(SizedMatrix{nq,nq,T})
∇qs   = rand(SizedMatrix{nq,nq,T})
∇us   = rand(SizedMatrix{nq,nu,T})
∇γs   = rand(SizedMatrix{nq,nγ,T})
∇bs   = rand(SizedMatrix{nq,nb,T})
∇q1s  = rand(SizedMatrix{nq,nq,T})


"""
	generate_dynamics_expressions(model::ContactDynamicsModel)
Generate fast dynamics methods using ModelingToolkit symbolic computing tools.
"""
function generate_dynamics_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b
	# Declare variables
	@variables q_1[1:nq]
	@variables q[1:nq]
	@variables u[1:nu]
	@variables γ[1:nγ]
	@variables b[1:nb]
	@variables q1[1:nq]
	@variables q̇[1:nq]
	@variables dt[1:1]

	# Symbolic computing
	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = ModelingToolkit.simplify.(L)

	dLq = ModelingToolkit.gradient(L, q, simplify=true)
	dLq̇ = ModelingToolkit.gradient(L, q̇, simplify=true)
	ddL = ModelingToolkit.sparsehessian(L, [q;q̇], simplify=true)
	ddL = SparseMatrixCSC{Expression,Int64}(ddL)
	ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

	# Mass Matrix
	M = M_func(model, q)
	M = ModelingToolkit.simplify.(M)

	# Control input Jacobian
	B = B_func(model, q)
	B = ModelingToolkit.simplify.(B)

	# Impact force Jacobian
	N = N_func(model, q)
	N = ModelingToolkit.simplify.(N)

	# Friction Force Jacobian
	P = _P_func(model, q)
	P = ModelingToolkit.simplify.(P)

	# Coriolis and Centrifugal forces Jacobians
	C = ddLq̇q * q̇ - dLq
	C = ModelingToolkit.simplify.(C)

	# Dynamics
	d = dynamics(model,dt[1],q_1,q,u,γ,b,q1)
	d = ModelingToolkit.simplify.(d)
	dz   = ModelingToolkit.jacobian(d, [q_1; q; u; γ; b; q1], simplify=true)
	dq_1 = ModelingToolkit.jacobian(d, q_1, simplify=true)
	dq   = ModelingToolkit.jacobian(d, q,   simplify=true)
	du   = ModelingToolkit.jacobian(d, q,   simplify=true)
	dγ   = ModelingToolkit.jacobian(d, γ,   simplify=true)
	db   = ModelingToolkit.jacobian(d, b,   simplify=true)
	dq1  = ModelingToolkit.jacobian(d, q1,  simplify=true)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = build_function(transpose([L]), q, q̇)[1] # need to transpose to get a line vector
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:N]    = build_function(N, q)[1]
	expr[:P]    = build_function(P, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	expr[:d]    = build_function(d,    dt, q_1, q, u, γ, b, q1)[1]
	expr[:dz]   = build_function(dz,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:dq_1] = build_function(dq_1, dt, q_1, q, u, γ, b, q1)[2]
	expr[:dq]   = build_function(dq,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:du]   = build_function(du,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:dγ]   = build_function(dγ,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:db]   = build_function(db,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:dq1]  = build_function(dq1,  dt, q_1, q, u, γ, b, q1)[2]
	return expr
end

# expr = generate_dynamics_expressions(model)
save_dynamics_expressions(expr, joinpath(@__DIR__, "quadruped_expr.jld2"), overwrite=true)
load_dynamics_expressions(joinpath(@__DIR__, "quadruped_expr.jld2"))
instantiate_dynamics!(model, joinpath(@__DIR__, "quadruped_expr.jld2"))





# L_fct, M_fct, B_fct, N_fct, P_fct, C_fct,
# 	d_fct, dz_fct, dq_1_fct, dq_fct, du_fct, dγ_fct, db_fct, dq1_fct = generate_dynamics(model)

dz_fct = generate_dynamics(model)
# dz_fct, dz2_fct, dz3_fct, dz4_fct = generate_dynamics(model)

z1  = [q1; q2; u2; γ2; b2; q3]
z1s = [q1s; q2s; u2s; γ2s; b2s; q3s]
nz = length(z1)
∇d1  = rand(SMatrix{nq,nz,T,nq*nz})
∇d1s = rand(SizedMatrix{nq,nz,T})

@btime dz_fct(model.dt, z1)
@btime dz_fct(model.dt, z1s)
ddd = dz_fct(∇d1s, model.dt, q1s, q2s, u2s, γ2s, b2s, q3s)
size(∇d1s)
sparse(∇d1s)
189/583
53*11
a = 10
a = 10
a = 10
a = 10

@btime ∇q1_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
@btime ∇q1_dynamics(model, model.dt, q1s, q2s, u2s, γ2s, b2s, q3s)

@btime dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
@btime dynamics(model, model.dt, q1s, q2s, u2s, γ2s, b2s, q3s)

@btime d_fct(model.dt, q1, q2, u2, γ2, b2, q3)
@btime d_fct(model.dt, q1s, q2s, u2s, γ2s, b2s, q3s)

@btime L_fct(q1s, q̇1s)
@btime M_fct(q1s)
@btime B_fct(q1s)
@btime N_fct(q1s)
@btime P_fct(q1s)
@btime C_fct(q1s, q̇1s)
@btime d_fct(model.dt, q1s, q2s, u2s, γ2s, b2s, q3s)
@btime dynamics(model, model.dt, q1s, q2s, u2s, γ2s, b2s, q3s)

@btime d_fct(model.dt, q1, q2, u2, γ2, b2, q3)
@btime dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)


a = 1; a += 2
a
b




#NOTE: if a new model is instantiated, re-run the lines below
@variables z_sym[1:quadruped.n]
l(z) = lagrangian(quadruped, view(z, 1:quadruped.nq), view(z, quadruped.nq .+ (1:quadruped.nq)))
_l = ModelingToolkit.simplify.(l(z_sym))
_dL = ModelingToolkit.gradient(_l, z_sym)
_dLq = view(_dL, 1:quadruped.nq)
_dLq̇ = view(_dL, quadruped.nq .+ (1:quadruped.nq))
_ddL = ModelingToolkit.sparsehessian(_l, z_sym)
_ddLq̇q = view(_ddL, quadruped.nq .+ (1:quadruped.nq), 1:quadruped.nq)

dL = eval(ModelingToolkit.build_function(_dL, z_sym)[1])
dLq = eval(ModelingToolkit.build_function(_dLq, z_sym)[1])
dLq̇ = eval(ModelingToolkit.build_function(_dLq̇, z_sym)[1])
ddLq̇q = eval(ModelingToolkit.build_function(_ddLq̇q, z_sym)[1])
ddL = eval(ModelingToolkit.build_function(_ddL, z_sym)[1])

function C_func(model::QuadrupedBasic11, q, q̇)
	ddLq̇q([q; q̇]) * q̇ - dLq(Vector([q; q̇]))
end





















"""
	generate_C_func(model::ContactDynamicsModel)
Generate a fast method for Coriolis and Centrifugal Dynamics terms computation using ModelingToolkit symbolic computing tools.
"""
function generate_C_func(model::ContactDynamicsModel)
	# Declare variables
	@variables q[1:model.dim.q]
	@variables q̇[1:model.dim.q]
	# Symbolic computing
	C = C_func(model, q, q̇)
	C = ModelingToolkit.simplify.(C)
	# Build function
	expr = build_function(C, q, q̇)[1]
	fct = eval(expr)
	return fct, expr
end

lagrangian_fct, lagrangian_expr = generate_lagrangian(model)
M_fct, M_expr = generate_M_func(model)
B_fct, B_expr = generate_B_func(model)
N_fct, N_expr = generate_N_func(model)
P_fct, P_expr = generate_P_func(model)
C_fct, C_expr = generate_C_func(model)

@time M_func(model, q1s)
@btime M_fct(M1s,q1s)


@btime B_func(model, q1s)
@btime B_fct(q1s)

SizedMatrix{nq,nq,T}(rand(nq,nq))



@btime lagrangian(model, q1, q̇1)



q1q̇1 = [q1;q̇1]
q1sq̇1s = [q1s;q̇1s]
@btime lagrangian_fct(q1, q̇1)







using ModelingToolkit, LinearAlgebra, SparseArrays

# Define the constants for the PDE
const α₂ = 1.0
const α₃ = 1.0
const β₁ = 1.0
const β₂ = 1.0
const β₃ = 1.0
const r₁ = 1.0
const r₂ = 1.0
const _DD = 100.0
const γ₁ = 0.1
const γ₂ = 0.1
const γ₃ = 0.1
const NN = 32
const X = reshape([i for i in 1:NN for j in 1:NN],NN,NN)
const Y = reshape([j for i in 1:NN for j in 1:NN],NN,NN)
const α₁ = 1.0.*(X.>=4*NN/5)

const Mx = Array(Tridiagonal([1.0 for i in 1:NN-1],[-2.0 for i in 1:NN],[1.0 for i in 1:NN-1]))
const My = copy(Mx)
Mx[2,1] = 2.0
Mx[end-1,end] = 2.0
My[1,2] = 2.0
My[end,end-1] = 2.0

# Define the discretized PDE as an ODE function
function f(u,p,t)
    A = u[:,:,1]
    B = u[:,:,2]
    C = u[:,:,3]
    MyA = My*A
    AMx = A*Mx
    DA = @. _DD*(MyA + AMx)
    dA = @. DA + α₁ - β₁*A - r₁*A*B + r₂*C
    dB = @. α₂ - β₂*B - r₁*A*B + r₂*C
    dC = @. α₃ - β₃*C + r₁*A*B - r₂*C
    cat(dA,dB,dC,dims=3)
end





# Define the initial condition as normal arrays
@variables u[1:NN,1:NN,1:3]
du = simplify.(f(u,nothing,0.0))
@show typeof(du)
@show typeof(u)
@show size(du)
@show size(u)

fastf = eval(ModelingToolkit.build_function(du,u,
            parallel=ModelingToolkit.MultithreadedForm())[2])



jac = ModelingToolkit.sparsejacobian(vec(du),vec(u))
fjac = eval(ModelingToolkit.build_function(jac,u,
            parallel=ModelingToolkit.MultithreadedForm())[2])

njac = eval(ModelingToolkit.build_function(jac,u)[2])

ModelingToolkit.build_function(jac,u,
            parallel=ModelingToolkit.MultithreadedForm())


using OrdinaryDiffEq
u0 = zeros(N,N,3)
MyA = zeros(N,N);
AMx = zeros(N,N);
DA = zeros(N,N);
prob = ODEProblem(f!,u0,(0.0,10.0))
fastprob = ODEProblem(ODEFunction((du,u,p,t)->fastf(du,u),
                                   jac = (du,u,p,t) -> fjac(du,u),
                                   jac_prototype = similar(jac,Float64)),
                                   u0,(0.0,10.0))


using BenchmarkTools
@btime solve(prob, TRBDF2()) # 33.073 s (895404 allocations: 23.87 GiB)
@btime solve(fastprob, TRBDF2()) # 209.670 ms (8208 allocations: 109.25 MiB)





















@time kinematics_1(model, q1; body = :torso, mode = :ee)
@btime kinematics_1(model, q1; body = :torso, mode = :ee)



function generate_kin1(model::QuadrupedModel)
	@variables q[1:model.dim.q]
	k1 = kinematics_1(model, q, body=:torso, mode=:ee)
	k1 = transpose(k1)
    k1 = ModelingToolkit.simplify.(k1)

	expr = build_function(k1, q)[1]
	fct = eval(expr)
    return fct, expr
end


mtk_kinematics_1, kin1_expr = generate_kin1(model)

k10 = kinematics_1(model, q1, body=:torso, mode=:ee)
@btime k11 = mtk_kinematics_1(q1)
@test norm(k10 .- k11') == 0.0

save("/home/simon/.julia/dev/Citron/proto_v5/dynamics_benchmark/kin1.jld2",
	"kin1_expr", kin1_expr, )

kin1_expr_loaded = load("/home/simon/.julia/dev/Citron/proto_v5/dynamics_benchmark/kin1.jld2",
	"kin1_expr", kin1_expr, )

@load "/home/simon/.julia/dev/Citron/proto_v5/dynamics_benchmark/kin1.jld2" kin1_expr

kin1_expr

@time loaded_mtk_kinematics_1 = eval(kin1_expr)



@btime dynamics_eval(model, model.dt, q1, q2, u2, γ2, b2, q3)
@ballocated dynamics_eval(model, model.dt, q1, q2, u2, γ2, b2, q3)










function dynamics_eval(model::ContactDynamicsModel, dt::T, qk_1::Vq_1, qk::Vq, uk::Vu,
	γk::Vγ, bk::Vb, qk1::Vq) where {T,Vq_1,Vq,Vu,Vγ,Vb,Vq1}

	v = Vector((qk1 - qk) / dt)
	joint_fric = model.joint_friction * v
	joint_fric[1:2] .= 0.0

	return ((1.0 / dt) * (M_func(model, qk_1) * (qk - qk_1)
	- M_func(model, qk) * (qk1 - qk))
	+ transpose(B_func(model, qk1)) * uk * 1.0
	+ transpose(N_func(model, qk1)) * γk * 1.0
	+ transpose(_P_func(model, qk1)) * bk * 1.0
	- dt * joint_fric
	- dt * (1.0 * C_func(model, qk1, (qk1 - qk) / dt)))

    return nothing
end


#
#
# function dynamics(model::QuadrupedModel15, dt::T, qk_1::SVq_1, qk::SVq,
# 	uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}
#
# 	v = Vector((qk1 - qk) / dt)
# 	joint_fric = model.joint_friction * v
# 	joint_fric[1:2] .= 0.0
#
# 	return ((1.0 / dt) * (M_func(model, qk_1) * (qk - qk_1)
# 	- M_func(model, qk) * (qk1 - qk))
# 	+ transpose(B_func(model, qk1)) * uk * 1.0
# 	+ transpose(N_func(model, qk1)) * γk * 1.0
# 	+ transpose(_P_func(model, qk1)) * bk * 1.0
# 	- dt * joint_fric
# 	- dt * (1.0 * C_func(model, qk1, (qk1 - qk) / dt)))
# end

function ∇q_1_dynamics(model::QuadrupedModel, dt::T, qk_1::SVq_1, qk::SVq,
	uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}
	dynx(x) = dynamics(model, dt, x, qk, uk, γk, bk, qk1)
	return ForwardDiff.jacobian(dynx, qk_1)
end

function ∇q_dynamics(model::QuadrupedModel, dt::T, qk_1::SVq_1, qk::SVq,
	uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}

	dynx(x) = dynamics(model, dt, qk_1, x, uk, γk, bk, qk1)
	return ForwardDiff.jacobian(dynx, qk)
end

function ∇u_dynamics(model::QuadrupedModel, dt::T, qk_1::SVq_1, qk::SVq,
	uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}

	dynx(x) = dynamics(model, dt, qk_1, qk, x, γk, bk, qk1)
	return ForwardDiff.jacobian(dynx, uk)
end

function ∇γ_dynamics(model::QuadrupedModel, dt::T, qk_1::SVq_1, qk::SVq,
	uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}

	dynx(x) = dynamics(model, dt, qk_1, qk, uk, x, bk, qk1)
	return ForwardDiff.jacobian(dynx, γk)
end

function ∇b_dynamics(model::QuadrupedModel, dt::T, qk_1::SVq_1, qk::SVq,
	uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}

	dynx(x) = dynamics(model, dt, qk_1, qk, uk, γk, x, qk1)
	return ForwardDiff.jacobian(dynx, bk)
end


function ∇q1_dynamics(model::QuadrupedModel, dt::T, qk_1::SVq_1, qk::SVq,
    uk::SVu, γk::SVγ, bk::SVb, qk1::SVq1) where {T,SVq_1,SVq,SVu,SVγ,SVb,SVq1}

	dynx(x) = dynamics(model, dt, qk_1, qk, uk, γk, bk, x)
	return ForwardDiff.jacobian(dynx, qk1)
end



# nq = 11
# nc = 4
# nγ = nc
# nb = 2nc
# nu = 8
#
# q1 = rand(SVector{nq,T})
# q2 = rand(SVector{nq,T})
# u2 = rand(SVector{nu,T})
# γ2 = rand(SVector{nγ,T})
# b2 = rand(SVector{nb,T})
# q3 = rand(SVector{nq,T})
#
# model

# _P_func(model, q3)
# @time dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
# @time ∇q_1_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
# @time ∇q_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
# @time ∇u_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
# @time ∇γ_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
# @time ∇b_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
# @time ∇q1_dynamics(model, model.dt, q1, q2, u2, γ2, b2, q3)
