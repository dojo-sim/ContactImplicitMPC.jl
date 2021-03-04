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

∇ys
∇q0s
∇q1s
∇u1s
∇γ1s
∇b1s
∇q2s



zs = rand(SizedVector{nz})
zs0 = deepcopy(zs)
θs = rand(SizedVector{nθ})
θs0 = deepcopy(θs)
ρs = 1e-3

residual(model, zs, θs, ρs)

zs1 = deepcopy(zs)
θs1 = deepcopy(θs)

zs1 == zs0
θs1 == θs0








function generate_residual_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b
	nz = model.dim.z
	nθ = model.dim.θ

	# Declare variables
	@variables dt[1:1]
	@variables z[1:nz]
	@variables θ[1:nθ]
	@variables ρ[1:1]
	# # Residual
	r = residual(model, dt[1], z, θ, ρ[1])
	# r = ModelingToolkit.simplify.(r)
	# rz = ModelingToolkit.sparsejacobian(r, z, simplify=true)
	# rθ = ModelingToolkit.sparsejacobian(r, θ, simplify=true)


	# Residual
	# r = [z[1], z[2], θ[1], ρ[1], dt[1]]
	r = ModelingToolkit.simplify.(r)
	rz = ModelingToolkit.sparsejacobian(r, z, simplify=true)
	rθ = ModelingToolkit.sparsejacobian(r, θ, simplify=true)


	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r,  dt, z, θ, ρ)[2]
	expr[:rz] = build_function(rz, dt, z, θ, ρ)[2]
	# expr[:rθ] = build_function(rθ, dt, z, θ, ρ)[2]
	return expr, r, rz, rθ
end


model = deepcopy(quadruped)
# instantiate_base!(model, joinpath(@__DIR__, ".expr/quadruped_base.jld2"))
# instantiate_dynamics!(model, joinpath(@__DIR__, ".expr/quadruped_dynamics.jld2"))
expr_res, r_, rz_, rθ_ = generate_residual_expressions(model)
rz_
rθ_
rθ_.nzval[1:10]
rθ_.nzval[11:20]
rθ_.nzval[21:30]
rθ_.nzval[31:39]
rθ_.colptr
rθ_.rowval







save_expressions(expr_res, joinpath(@__DIR__, ".expr/quadruped_residual_dummy.jld2"), overwrite=true)
instantiate_residual!(model, joinpath(@__DIR__, ".expr/quadruped_residual_dummy.jld2"))


aaa = [1 0 0 0 2;
	   3 4 5 0 1]
bbb = sparse(aaa)



zs = rand(SizedVector{nz})
θs = rand(SizedVector{nθ})
ρs = 1e-3
rs = rand(SizedVector{nz})
∇zs = rand(nz,nz)
∇θs = rand(nz,nθ)
@btime r_fast!(model, rs, zs, θs, ρs)
@btime rz_fast!(model, ∇zs, zs, θs, ρs)
@btime rθ_fast!(model, ∇θs, zs, θs, ρs)

∇zs
sparse(∇zs)
215/nz^2
sparse(∇θs)
162/(nz*nθ)
