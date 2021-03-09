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
	h =   θ[off .+ (1:1)]
	return q0, q1, u1, w1, h
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
	ψ =  z[off .+ (1:nc)]
	off += nc
	η =  z[off .+ (1:nb)]
	off += nb
	s1 = z[off .+ (1:nc)]
	off += nc
	s2 = z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, ψ, η, s1, s2
end

function z_initialize!(z, model::ContactDynamicsModel, q1)
	nq = model.dim.q

    z .= 1.0
    z[1:nq] = q1
end

function θ_initialize!(θ, model::ContactDynamicsModel, q0, q1, u, w, h)
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
    θ[off .+ (1:1)] .= h
end

function num_var(model::ContactDynamicsModel)
    model.dim.q + model.dim.c + model.dim.b + model.dim.c + model.dim.b + 2 * model.dim.c
end

function num_data(model::ContactDynamicsModel)
    model.dim.q + model.dim.q + model.dim.u + model.dim.w + 1
end

inequality_indices(model::ContactDynamicsModel) = collect(model.dim.q .+ (1:(num_var(model) - model.dim.q)))

function dynamics(model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)

	v1 = (q2 - q1) / h[1]

	return (1.0 / h[1] *
		  (M_fast(model, q0) * (q1 - q0)
		- M_fast(model, q1) * (q2 - q1))
		+ transpose(B_fast(model, q2)) * u1
		+ transpose(A_fast(model, q2)) * w1
		+ transpose(N_fast(model, q2)) * γ1
		+ transpose(P_fast(model, q2)) * b1
		- h[1] * C_fast(model, q2, v1)
		- h[1] * model.joint_friction .* v1)
end

function residual(model::ContactDynamicsModel, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b

	q0, q1, u1, w1, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ, η, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)
	vT = (P_fast(model, q2) * q2 - P_fast(model, q1) * q1) / h[1]

	[dynamics(model, h, q0, q1, u1, w1, γ1, b1, q2);
	 s1 - ϕ;
	 γ1 .* s1 .- κ;
	 vT + vcat([ψi .* ones(Int(nb / nc)) for ψi in ψ]...) - η;
	 s2 .- (model.μ_world * γ1
	   .- transpose(sum(reshape(b1, (Int(nb/nc), nc)), dims = 1))[:, 1]);
	 ψ .* s2 .- κ;
	 b1 .* η .- κ]
end

mutable struct BaseMethods
	L::Any
	M::Any
	B::Any
	A::Any
	N::Any
	P::Any
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
	dq0::Any
	dq1::Any
	du1::Any
	dw1::Any
	dγ1::Any
	db1::Any
	dq2::Any
end

function DynamicsMethods()
	function f()
		error("Not Implemented: use instantiate_dynamics!")
		return nothing
	end
	return DynamicsMethods(fill(f, 9)...)
end

mutable struct ResidualMethods
	r::Any
	rz::Any
	rθ::Any
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
end

```
	get_model(name::String)
	Helper function that provides a model where fast functions have been instantiated.
```
function get_model(name::String)
	path = joinpath(@__DIR__, name)
	include(joinpath(path, "model.jl"))
	if name == "particle"
		model = particle
	elseif name == "quadruped"
		model = quadruped
	else
		error("Unknown model name")
	end
	instantiate_base!(model, joinpath(path, "base.jld2"))
	instantiate_dynamics!(model, joinpath(path, "dynamics.jld2"))
	instantiate_residual!(model, joinpath(path, "residual.jld2"))
	@load joinpath(path, "sparse_jacobians.jld2") rz_sp rθ_sp
	model.spa.rz_sp = rz_sp
	return model
end
