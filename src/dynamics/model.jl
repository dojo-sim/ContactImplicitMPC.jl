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

function num_var(model)
    model.dim.q + model.dim.c + model.dim.b + model.dim.c + model.dim.b + 2 * model.dim.c
end

function num_data(model)
    model.dim.q + model.dim.q + model.dim.u + model.dim.w + 1
end

function dynamics(model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)

	v1 = (q2 - q1) / h[1]
	joint_friction = [zeros(2); model.μ_joint * v1[3:end]] # TODO: only good for planar systems

	return (1.0 / h[1] *
		  (M_func(model, q0) * (q1 - q0)
		- M_func(model, q1) * (q2 - q1))
		+ transpose(B_func(model, q2)) * u1
		+ transpose(N_func(model, q2)) * γ1
		+ transpose(P_func(model, q2)) * b1
		- h[1] * C_func(model, q2, v1)
		- h[1] * joint_friction)
end

function residual(model::ContactDynamicsModel, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b

	q0, q1, u1, w1, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ, η, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)
	vT = (P_func(model, q2) * q2 - P_func(model, q1) * q1) / h[1]

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
	N::Any
	P::Any
	C::Any
end

function BaseMethods()
	function f()
		error("Not Implemented: use instantiate_base!")
		return nothing
	end
	return BaseMethods(fill(f, 6)...)
end

mutable struct DynamicsMethods
	d::Any
	dy::Any
	dq0::Any
	dq1::Any
	du1::Any
	dγ1::Any
	db1::Any
	dq2::Any
end

function DynamicsMethods()
	function f()
		error("Not Implemented: use instantiate_dynamics!")
		return nothing
	end
	return DynamicsMethods(fill(f, 8)...)
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
