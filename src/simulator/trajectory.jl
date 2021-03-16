struct ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}
	H::Int
	h::T
	q::Vector{SizedArray{Tuple{nq},T,1,1}}         # trajectory of q's   length=H+2
	u::Vector{SizedArray{Tuple{nu},T,1,1}}         # trajectory of u's   length=H
	w::Vector{SizedArray{Tuple{nw},T,1,1}}         # trajectory of w's   length=H
	γ::Vector{SizedArray{Tuple{nc},T,1,1}}         # trajectory of γ's   length=H
	b::Vector{SizedArray{Tuple{nb},T,1,1}}         # trajectory of b's   length=H
	z::Vector{SizedArray{Tuple{nz},T,1,1}}         # trajectory of z's   length=H
	θ::Vector{SizedArray{Tuple{nθ},T,1,1}}         # trajectory of θ's   length=H
	κ::Vector{T}
end

function contact_trajectory(H::Int, h::T, model::ContactDynamicsModel) where {T}
	dim = model.dim
	nq = dim.q
    nu = dim.u
    nw = dim.w
    nc = dim.c
    nb = dim.b
	nz = num_var(dim)
	nθ = num_data(dim)

	q = [zeros(SizedVector{nq,T}) for k=1:H+2]
	u = [zeros(SizedVector{nu,T}) for k=1:H]
	w = [zeros(SizedVector{nw,T}) for k=1:H]
	γ = [zeros(SizedVector{nc,T}) for k=1:H]
	b = [zeros(SizedVector{nb,T}) for k=1:H]
	z = [zeros(SizedVector{nz,T}) for k=1:H]
	θ = [zeros(SizedVector{nθ,T}) for k=1:H]
	κ = [0.0]

	return ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}(H,h,q,u,w,γ,b,z,θ,κ)
end

struct ContactDerivTraj{S,nq,nu,nc,nb}
	dq2dq0::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
    dq2dq1::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
    dq2du::Vector{SizedArray{Tuple{nq,nu},S,2,2}}
    dγdq0::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
    dγdq1::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
    dγdu::Vector{SizedArray{Tuple{nc,nu},S,2,2}}
    dbdq0::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
    dbdq1::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
    dbdu::Vector{SizedArray{Tuple{nb,nu},S,2,2}}
end

function contact_derivative_trajectory(H::Int, model::ContactDynamicsModel)
	nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = model.dim.b

	dq2dq0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:H]
	dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:H]
	dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:H]
	dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:H]
	dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:H]
	dγdu = [SizedMatrix{nc,nu}(zeros(nc, nu)) for t = 1:H]
	dbdq0 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:H]
	dbdq1 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:H]
	dbdu = [SizedMatrix{nb,nu}(zeros(nb, nu)) for t = 1:H]

	ContactDerivTraj(
		dq2dq0,
		dq2dq1,
		dq2du,
		dγdq0,
		dγdq1,
		dγdu,
		dbdq0,
		dbdq1,
		dbdu)
end
