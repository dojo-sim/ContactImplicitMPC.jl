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
	y::Int         # whole dim
end

function Dimensions12(q::Int,u::Int,c::Int,b::Int)
    γ = c
	n = 2q
	θ = 2q + u
	w = γ + b + q
    y = 3q + u + γ + b
    return Dimensions12(q,u,γ,b,c,n,θ,w,y)
end

mutable struct Indices13{Vq_1,Vq,Vu,Vγ,Vb,Vq1,Vθ,Vw}
    q0::Vq_1      # qk-1
    q1::Vq        # qk
    u1::Vu        # uk
    γ1::Vγ        # γk
    b1::Vb        # bk
    q2::Vq1       # qk+1
	θ::Vθ         # θ = [qk-1; qk; u_k]
	w::Vw         # w = [γk; bk; qk+1]
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

mutable struct ResidualMethods12
	r::Any
	rz::Any
	rθ::Any
end

function ResidualMethods12()
	function f()
		error("Not Implemented: use instantiate_residual!")
		return nothing
	end
	return ResidualMethods12(fill(f, 3)...)
end
