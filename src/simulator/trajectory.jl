mutable struct ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}
	H::Int                                         # horizon length
	h::T                                           # time step length
	q::Vector{SizedArray{Tuple{nq},T,1,1}}         # trajectory of q's   length=H+2
	u::Vector{SizedArray{Tuple{nu},T,1,1}}         # trajectory of u's   length=H
	w::Vector{SizedArray{Tuple{nw},T,1,1}}         # trajectory of w's   length=H
	γ::Vector{SizedArray{Tuple{nc},T,1,1}}         # trajectory of γ's   length=H
	b::Vector{SizedArray{Tuple{nb},T,1,1}}         # trajectory of b's   length=H
	z::Vector{SizedArray{Tuple{nz},T,1,1}}         # trajectory of z's   length=H
	θ::Vector{SizedArray{Tuple{nθ},T,1,1}}         # trajectory of θ's   length=H
	κ::T
end

function ContactTraj(H::Int, h::T, model::ContactDynamicsModel) where {T}
	return ContactTraj(H, h, model.dim)
end

function ContactTraj(H::Int, h::T, dim::Dimensions) where {T}
	nz = num_var(dim)
	nθ = num_data(dim)
	return ContactTraj(H, h, dim.q, dim.u, dim.w, dim.c, dim.b, nz, nθ)
end

function ContactTraj(H::Int, h::T, nq::Int, nu::Int, nw::Int, nc::Int, nb::Int,
	nz::Int, nθ::Int) where{T}
	q = [zeros(SizedVector{nq,T}) for k=1:H+2]
	u = [zeros(SizedVector{nu,T}) for k=1:H]
	w = [zeros(SizedVector{nw,T}) for k=1:H]
	γ = [zeros(SizedVector{nc,T}) for k=1:H]
	b = [zeros(SizedVector{nb,T}) for k=1:H]
	z = [zeros(SizedVector{nz,T}) for k=1:H]
	θ = [zeros(SizedVector{nθ,T}) for k=1:H]
	κ = 0.0
	return ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}(H,h,q,u,w,γ,b,z,θ,κ)
end
