abstract type ContactDynamicsModel
end

abstract type QuadrupedModel <: ContactDynamicsModel
end

mutable struct Dimensions12
    q::Int         # configuration dim
    u::Int         # control dim
	γ::Int         # normal components
	b::Int         # linear friction components
    b̄::Int         # minimal friction components
	c::Int         # number of contact points
	n::Int         # state dim
	θ::Int         # data dim
	w::Int         # sol dim
	y::Int         # dynamics variables dim
	z::Int         # z dim
end

function Dimensions12(q::Int,u::Int,c::Int,b::Int)
    γ = c
	b̄ = Int(b/2)
	n = 2q
	θ = 2q + u
	w = γ + b + q
    y = 3q + u + γ + b
	z = q + 4γ + 2b
    return Dimensions12(q,u,γ,b,b̄,c,n,θ,w,y,z)
end

mutable struct Indices13{Vq_1,Vq,Vu,Vγ,Vb,Vq1,Vθ,Vw}
    q0::Vq_1      # q1-1
    q1::Vq        # q1
    u1::Vu        # uk
    γ1::Vγ        # γk
    b1::Vb        # bk
    q2::Vq1       # q1+1
	θ::Vθ         # θ = [q1-1; q1; u_k]
	w::Vw         # w = [γk; bk; q1+1]
end

function Indices13(q::Int,u::Int,c::Int,b::Int)
	off = 0
	iq0 = SizedVector{q,Int}(off .+ (1:q)); off += q;
	iq1 = SizedVector{q,Int}(off .+ (1:q)); off += q
	iu1 = SizedVector{u,Int}(off .+ (1:u)); off += u
	iγ1 = SizedVector{c,Int}(off .+ (1:c)); off += c
	ib1 = SizedVector{b,Int}(off .+ (1:b)); off += b
	iq2 = SizedVector{q,Int}(off .+ (1:q)); off += q
	iθ  = [iq0; iq1; iu1]
	iw  = [iγ1; ib1; iq2]
    return Indices13(iq0,iq1,iu1,iγ1,ib1,iq2,iθ,iw)
end

function unpack_θ(model::ContactDynamicsModel, θ::AbstractVector{T}) where T
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b

	# Parameters
	off = 0
	q0 = θ[off .+ (1:nq)]
	off += nq
	q1 = θ[off .+ (1:nq)]
	off += nq
	u1 =  θ[off .+ (1:nu)]
	return q0, q1, u1
end

function unpack_z(model::ContactDynamicsModel, z::AbstractVector{T}) where {T}
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b

	# system variables
	off = 0
	q2 =  z[off .+ (1:nq)]
	off += nq
	γ1 =  z[off .+ (1:nγ)]
	off += nγ
	b1 =  z[off .+ (1:nb)]
	off += nb
	ψ =  z[off .+ (1:nγ)]
	off += nγ
	η =  z[off .+ (1:nb)]
	off += nb
	s1 = z[off .+ (1:nγ)]
	off += nγ
	s2 = z[off .+ (1:nγ)]
	off += nγ
	return q2, γ1, b1, ψ, η, s1, s2
end

mutable struct BaseMethods12
	L::Any
	M::Any
	B::Any
	N::Any
	P::Any
	C::Any
end

function BaseMethods12()
	function f()
		error("Not Implemented: use instantiate_base!")
		return nothing
	end
	return BaseMethods12(fill(f, 6)...)
end

mutable struct DynamicsMethods13
	d::Any
	dy::Any
	dq0::Any
	dq1::Any
	du1::Any
	dγ1::Any
	db1::Any
	dq2::Any
end

function DynamicsMethods13()
	function f()
		error("Not Implemented: use instantiate_dynamics!")
		return nothing
	end
	return DynamicsMethods13(fill(f, 8)...)
end

mutable struct ResidualMethods14
	r::Any
	rz::Any
	rθ::Any
end

function ResidualMethods14()
	function f()
		error("Not Implemented: use instantiate_residual!")
		return nothing
	end
	return ResidualMethods14(fill(f, 3)...)
end
