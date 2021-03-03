abstract type ContactDynamicsModel
end

abstract type QuadrupedModel <: ContactDynamicsModel
end

mutable struct Dimensions12
    q::Int         # configuration dim
    u::Int         # control dim
	γ::Int         # normal components
    b::Int         # friction components
	c::Int         # number of contact points
	n::Int         # state dim
	θ::Int         # data dim
	w::Int         # sol dim
	z::Int         # whole dim
end

function Dimensions12(q::Int,u::Int,c::Int,b::Int)
    γ = c
	n = 2q
	θ = 2q + u
	w = γ + b + q
    z = 3q + u + γ + b
    return Dimensions12(q,u,γ,b,c,n,θ,w,z)
end

mutable struct Indices13{Vq_1,Vq,Vu,Vγ,Vb,Vq1,Vθ,Vw}
    q_1::Vq_1     # qk-1
    q::Vq         # qk
    u::Vu         # uk
    γ::Vγ         # γk
    b::Vb         # bk
    q1::Vq1       # qk+1
	θ::Vθ         # θ = [qk-1; qk; u_k]
	w::Vw         # w = [γk; bk; qk+1]
end

function Indices13(q::Int,u::Int,c::Int,b::Int)
	off = 0
	iq_1 = SizedVector{q,Int}(off .+ (1:q)); off += q;
	iq   = SizedVector{q,Int}(off .+ (1:q)); off += q
	iu   = SizedVector{u,Int}(off .+ (1:u)); off += u
	iγ   = SizedVector{c,Int}(off .+ (1:c)); off += c
	ib   = SizedVector{b,Int}(off .+ (1:b)); off += b
	iq1  = SizedVector{q,Int}(off .+ (1:q)); off += q
	iθ   = [iq_1; iq; iu]
	iw   = [iγ; ib; iq1]
    return Indices13(iq_1,iq,iu,iγ,ib,iq1,iθ,iw)
end

mutable struct DynamicsMethods11
	L::Any
	M::Any
	B::Any
	N::Any
	P::Any
	C::Any
	d::Any
	dz::Any
	dq_1::Any
	dq::Any
	du::Any
	dγ::Any
	db::Any
	dq1::Any
end

function DynamicsMethods11()
	function f()
		return nothing
	end
	return DynamicsMethods11(fill(f, 14)...)
end
