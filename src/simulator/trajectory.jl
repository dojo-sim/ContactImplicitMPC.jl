mutable struct ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}
	H::Int                               # horizon length
	q::Vector{SizedVector{nq,T}}         # trajectory of q's   length=H+2
	u::Vector{SizedVector{nu,T}}         # trajectory of u's   length=H
	w::Vector{SizedVector{nw,T}}         # trajectory of w's   length=H
	γ::Vector{SizedVector{nc,T}}         # trajectory of γ's   length=H
	b::Vector{SizedVector{nb,T}}         # trajectory of b's   length=H
	z::Vector{SizedVector{nz,T}}         # trajectory of z's   length=H
	θ::Vector{SizedVector{nθ,T}}         # trajectory of θ's   length=H
	κ::T
end

function ContactTraj(H::Int, model::ContactDynamicsModel; T::DataType=Float64)
	nz = num_var(model)
	nθ = num_data(model)
	return ContactTraj(H, model.dim.q, model.dim.u, model.dim.w, model.dim.c, model.dim.b, nz, nθ, T=T)
end

function ContactTraj(H::Int, nq::Int, nu::Int, nw::Int, nc::Int, nb::Int,
	nz::Int, nθ::Int; T::DataType=Float64)

	q = [zeros(SizedVector{nq,T}) for k=1:H+2]
	u = [zeros(SizedVector{nu,T}) for k=1:H]
	w = [zeros(SizedVector{nw,T}) for k=1:H]
	γ = [zeros(SizedVector{nc,T}) for k=1:H]
	b = [zeros(SizedVector{nb,T}) for k=1:H]
	z = [zeros(SizedVector{nz,T}) for k=1:H]
	θ = [zeros(SizedVector{nθ,T}) for k=1:H]
	κ = 0.0
	return ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}(H,q,u,w,γ,b,z,θ,κ)
end
