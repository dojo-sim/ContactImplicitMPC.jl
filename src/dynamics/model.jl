abstract type ContactDynamicsModel end

struct Dimensions
    q::Int         # configuration
    u::Int         # control
	w::Int         # disturbance
	c::Int         # contact points
	b::Int         # linear friction components
end

function unpack_θ(model::ContactDynamicsModel, θ)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b

	# Parameters
	off = 0
	q0 = θ[off .+ (1:nq)]
	off += nq
	q1 = θ[off .+ (1:nq)]
	off += nq
	u1 =  θ[off .+ (1:nu)]
	off += nu
	w1 =  θ[off .+ (1:nw)]
	off += nw
	μ =   θ[off .+ (1:1)]
	off += 1
	h =   θ[off .+ (1:1)]
	return q0, q1, u1, w1, μ, h
end

function pack_θ(model::ContactDynamicsModel, q0, q1, u1, w1, μ, h)
	return [q0; q1; u1; w1; μ; h]
end

function unpack_z(model::ContactDynamicsModel, z)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b

	# system variables
	off = 0
	q2 =  z[off .+ (1:nq)]
	off += nq
	γ1 =  z[off .+ (1:nc)]
	off += nc
	b1 =  z[off .+ (1:nb)]
	off += nb
	ψ1 =  z[off .+ (1:nc)]
	off += nc
	η1 =  z[off .+ (1:nb)]
	off += nb
	s1 = z[off .+ (1:nc)]
	off += nc
	s2 = z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, ψ1, η1, s1, s2
end

function pack_z(model::ContactDynamicsModel, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_fast(model, q2)
	s2 = model.μ_world * γ1 .- E_func(model) * b1
	return [q2; γ1; b1; ψ1; η1; s1; s2]
end

function z_initialize!(z, model::ContactDynamicsModel, q1)
	nq = model.dim.q
    z .= 1.0
    z[1:nq] = q1
end

function θ_initialize!(θ, model::ContactDynamicsModel, q0, q1, u, w, μ, h)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	off = 0
    θ[1:nq] = q0
	off += nq
    θ[off .+ (1:nq)] = q1
	off += nq
    θ[off .+ (1:nu)] = u
	off += nu
    θ[off .+ (1:nw)] = w
	off += nw
	θ[off .+ (1:1)] .= μ
	off += 1
    θ[off .+ (1:1)] .= h
end

function num_var(model::ContactDynamicsModel)
	num_var(model.dim)
end

function num_data(model::ContactDynamicsModel)
	num_data(model.dim)
end

function num_var(dim::Dimensions)
    dim.q + dim.c + dim.b + dim.c + dim.b + 2 * dim.c
end

function num_data(dim::Dimensions)
    dim.q + dim.q + dim.u + dim.w + 1 + 1
end

inequality_indices(model::ContactDynamicsModel) = collect(model.dim.q .+ (1:(num_var(model) - model.dim.q)))

function E_func(model::ContactDynamicsModel)
	SMatrix{model.dim.c, model.dim.b}(kron(Diagonal(ones(model.dim.c)), ones(1, Int(model.dim.b / model.dim.c))))
end

# https://github.com/HarvardAgileRoboticsLab/drake/blob/75b260c9eb250d08ffbbf3daa80758e4fe558d7f/drake/matlab/solvers/trajectoryOptimization/VariationalTrajectoryOptimization.m
function lagrangian_derivatives(model::ContactDynamicsModel, q, v)
	D1L = -1.0 * C_fast(model, q, v)
    D2L = M_fast(model, q) * v
	return D1L, D2L
end

function dynamics(model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)

	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	ne = dim(model.env)
	k = kinematics(model, q2)
	λ1 = vcat([transpose(rotation(model.env, k[(i-1) * ne .+ (1:ne)])) * [friction_mapping(model.env) * b1[(i-1) * nf .+ (1:nf)]; γ1[i]] for i = 1:nc]...) # TODO: make efficient

	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1
		- h[1] * model.joint_friction .* vm2)
end

function dynamics(model::ContactDynamicsModel, h, q0, q1, u1, w1, λ1, q2)

	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1
		- h[1] * model.joint_friction .* vm2)
end

function contact_forces(model::ContactDynamicsModel, γ1, b1, q2)
	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	ne = dim(model.env)
	k = kinematics(model, q2)
	λ1 = vcat([transpose(rotation(model.env, k[(i-1) * ne .+ (1:ne)])) * [friction_mapping(model.env) * b1[(i-1) * nf .+ (1:nf)]; γ1[i]] for i = 1:nc]...) # TODO: make efficient
end

function velocity_stack(model::ContactDynamicsModel, q1, q2, h)
	nc = model.dim.c
	np = dim(model.env)
	k = kinematics(model, q2)
	v = J_fast(model, q2) * (q2 - q1) / h[1]
	v_surf = [rotation(model.env, k[(i-1) * np .+ (1:np)]) * v[(i-1) * np .+ (1:np)] for i = 1:nc]
	vT_stack = vcat([[v_surf[i][1]; -v_surf[i][1]] for i = 1:nc]...)
end

function residual(model::ContactDynamicsModel, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	np = dim(model.env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)

	λ1 = contact_forces(model, γ1, b1, q2)
	vT_stack = velocity_stack(model, q1, q2, h)
	ψ_stack = transpose(E_func(model)) * ψ1

	[model.dyn.d(h, q0, q1, u1, w1, λ1, q2);
	 vT_stack + ψ_stack - η1;
	 s1 - ϕ;
	 s2 .- (μ[1] * γ1 .- E_func(model) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

function res_con(model::ContactDynamicsModel, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	np = dim(model.env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)

	[s1 - ϕ_func(model, q2);
	 s2 .- (μ[1] * γ1 .- E_func(model) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

function linearization_var_index(model::ContactDynamicsModel)
	nq = model.dim.q
	nb = model.dim.b
	nc = model.dim.c
	ny = 2nc + nb
	# x = [q2]
	# y1 = [γ1, b1, ψ1]
	# y2 = [s1, η1, s2]
	off = 0
	ix  = off .+ Vector(1:nq); off += nq
	iy1 = off .+ Vector(1:ny); off += ny
	iy2 = off .+ [Vector(nb .+ (1:nc)); Vector(1:nb); Vector(nb+nc .+ (1:nc))]; off += ny
	return ix, iy1, iy2
end

function linearization_term_index(model::ContactDynamicsModel)
	nq = model.dim.q
	nb = model.dim.b
	nc = model.dim.c
	ny = 2nc + nb
	# dyn = [dyn]
	# rst = [s1  - ..., ≡ ialt
	#        ... - η1 ,
	#        s2  - ...,]
	# bil = [γ1 .* s1 .- κ;
	#        b1 .* η1 .- κ;
	#        ψ1 .* s2 .- κ]
	idyn  = Vector(1:nq)
	irst2 = Vector(nq .+ (1:nb))
	irst1 = Vector(nq + nb .+ (1:nc))
	irst3 = Vector(nq + nb + nc .+ (1:nc))
	irst  = [irst1; irst2; irst3]
	ialt  = irst1
	ibil1 = Vector(nq + nb + 2nc .+ (1:nc))
	ibil2 = Vector(nq + nb + 3nc .+ (1:nb))
	ibil3 = Vector(nq + nb + 3nc + nb .+ (1:nc))
	ibil  = [ibil1; ibil2; ibil3]
	return idyn, irst, ibil, ialt
end

function get_bilinear_indices(model::ContactDynamicsModel)
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b

	terms = [SVector{nc,Int}(nq + nb + 2nc .+ (1:nc)),
			 SVector{nb,Int}(nq + nb + 3nc .+ (1:nb)),
			 SVector{nc,Int}(nq + nb + 3nc + nb .+ (1:nc))]

	vars = [[SVector{nc,Int}(nq .+ (1:nc)),           SVector{nc}(nq + 2nc + 2nb .+ (1:nc))], # γ1, s1
			[SVector{nb,Int}(nq + nc .+ (1:nb)),      SVector{nb}(nq + 2nc +  nb .+ (1:nb))], # b1, η
			[SVector{nc,Int}(nq + nc + nb .+ (1:nc)), SVector{nc}(nq + 3nc + 2nb .+ (1:nc))], # ψ, s2
			]
	return terms, vars
end

mutable struct BaseMethods
	L::Any
	ϕ::Any
	M::Any
	B::Any
	A::Any
	J::Any
	C::Any
end

function BaseMethods()
	function f()
		error("Not Implemented: use instantiate_base!")
		return nothing
	end
	return BaseMethods(fill(f, 7)...)
end

mutable struct DynamicsMethods
	d::Any
	dy::Any
	dθ::Any
	dq0::Any
	dq1::Any
	du1::Any
	dw1::Any
	dγ1::Any
	db1::Any
	dλ1::Any
	dq2::Any
end

function DynamicsMethods()
	function f()
		error("Not Implemented: use instantiate_dynamics!")
		return nothing
	end
	return DynamicsMethods(fill(f, 11)...)
end

function ResidualMethods()
	function f()
		error("Not Implemented: use instantiate_residual!")
		return nothing
	end
	return ResidualMethods(fill(f, 3)...)
end

mutable struct SparseStructure
	rz_sp::Any
	rθ_sp::Any
end

mutable struct ContactMethods
	cf
	dcf
	vs
	vsq2
	vsq1h
	mdvs
	mdψη
	rc
	rcz
	rcθ
end

function ContactMethods()
	function f()
		error("Not Implemented: use instantiate_contact_methods!")
		return nothing
	end
	return ContactMethods(fill(f, 10)...)
end

"""
	get_model(name::String, surf::String)
	Helper function that provides a model where fast functions have been instantiated.
"""
function get_model(name::String; model_name::String = name, surf::String = "flat")
	#TODO: assert model exists
	path = joinpath(@__DIR__, name)
	# include(joinpath(path, "model.jl"))
	model = eval(Symbol(model_name * (surf != "flat" ? "_" * surf : "")))
	instantiate_base!(model, joinpath(path, "dynamics", "base.jld2"))
	instantiate_dynamics!(model, joinpath(path, "dynamics", "dynamics.jld2"))
	instantiate_residual!(model, joinpath(path, surf, "residual.jld2"))
	instantiate_linearized!(model, joinpath(path, surf, "linearized.jld2"))
	@load joinpath(path, surf, "sparse_jacobians.jld2") rz_sp rθ_sp
	model.spa.rz_sp = rz_sp
	model.spa.rθ_sp = rθ_sp
	return model
end

function get_gait(name::String, gait::String)
	#TODO: assert model exists
	path = joinpath(@__DIR__, name)
	gait_path = joinpath(path, "gaits/" * gait * ".jld2")

	res = JLD2.jldopen(gait_path) # z̄ x̄ ū h̄ q u γ b

	return res["q"], res["u"], res["γ"], res["b"], mean(res["h̄"])
end

function get_trajectory(name::String, gait::String; model_name = name, load_type::Symbol=:split_traj,
	model::ContactDynamicsModel=eval(Symbol(model_name)))
	#TODO: assert model exists
	path = joinpath(@__DIR__, name)
	gait_path = joinpath(path, "gaits/" * gait * ".jld2")

	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b
	res = JLD2.jldopen(gait_path)

	if load_type == :split_traj
		q, u, γ, b, h = res["q"], res["u"], res["γ"], res["b"], mean(res["h̄"])
		ū = res["ū"]
		ψ = [ut[nu + nc + nb .+ (1:nc)] for ut in ū]
		η = [ut[nu + nc + nb + nc .+ (1:nb)] for ut in ū]

		T = length(u)

		traj = contact_trajectory(T, h, model)
		traj.q .= deepcopy(q)
		traj.u .= deepcopy(u)
		traj.γ .= deepcopy(γ)
		traj.b .= deepcopy(b)
		traj.z .= [pack_z(model, q[t+2], γ[t], b[t], ψ[t], η[t]) for t = 1:T]
		traj.θ .= [pack_θ(model, q[t], q[t+1], u[t], zeros(nw), model.μ_world, h) for t = 1:T]
	elseif load_type == :joint_traj
		traj = res["traj"]
	end

	return traj
end

function model_name(model::ContactDynamicsModel)
    name = Symbol(string(typeof(model).name)[10:end-1])
    return name
end
