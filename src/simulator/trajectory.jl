struct ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}
	H::Int
	h::T
	κ::Vector{T}
	q::Vector{SizedArray{Tuple{nq},T,1,1}}         # trajectory of q's   length=H+2
	u::Vector{SizedArray{Tuple{nu},T,1,1}}         # trajectory of u's   length=H
	w::Vector{SizedArray{Tuple{nw},T,1,1}}         # trajectory of w's   length=H
	γ::Vector{SizedArray{Tuple{nc},T,1,1}}         # trajectory of γ's   length=H
	b::Vector{SizedArray{Tuple{nb},T,1,1}}         # trajectory of b's   length=H
	z::Vector{SizedArray{Tuple{nz},T,1,1}}         # trajectory of z's   length=H
	θ::Vector{SizedArray{Tuple{nθ},T,1,1}}         # trajectory of θ's   length=H
	iq0::SizedArray{Tuple{nq},Int,1,1,Vector{Int}}
	iq1::SizedArray{Tuple{nq},Int,1,1,Vector{Int}}
	iu1::SizedArray{Tuple{nu},Int,1,1,Vector{Int}}
	iw1::SizedArray{Tuple{nw},Int,1,1,Vector{Int}}
	iq2::SizedArray{Tuple{nq},Int,1,1,Vector{Int}}
	iγ1::SizedArray{Tuple{nc},Int,1,1,Vector{Int}}
	ib1::SizedArray{Tuple{nb},Int,1,1,Vector{Int}}
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

	off = 0
	iq0 = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q0
	iq1 = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q1
	iu1 = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iw1 = SizedVector{nw}(off .+ (1:nw)); off += nw # index of the disturbance w1
	off = 0
	iq2 = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iγ1 = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib1 = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
	return ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}(H,h,κ,q,u,w,γ,b,z,θ,iq0,iq1,iu1,iw1,iq2,iγ1,ib1)
end

struct ContactDerivTraj{S,nq,nu,nc,nb}
	dq2dq0::Vector{SizedArray{Tuple{nq,nq},S,2,2,Array{S,2}}}
    dq2dq1::Vector{SizedArray{Tuple{nq,nq},S,2,2,Array{S,2}}}
    dq2du::Vector{SizedArray{Tuple{nq,nu},S,2,2,Array{S,2}}}
    dγdq0::Vector{SizedArray{Tuple{nc,nq},S,2,2,Array{S,2}}}
    dγdq1::Vector{SizedArray{Tuple{nc,nq},S,2,2,Array{S,2}}}
    dγdu::Vector{SizedArray{Tuple{nc,nu},S,2,2,Array{S,2}}}
    dbdq0::Vector{SizedArray{Tuple{nb,nq},S,2,2,Array{S,2}}}
    dbdq1::Vector{SizedArray{Tuple{nb,nq},S,2,2,Array{S,2}}}
    dbdu::Vector{SizedArray{Tuple{nb,nu},S,2,2,Array{S,2}}}
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

function update_z!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}, t::Int) where {T,nq,nu,nw,nc,nb,nz,nθ}
	if t in (1:traj.H)
		traj.z[t][traj.iq2] = traj.q[t+2]
		traj.z[t][traj.iγ1] = traj.γ[t]
		traj.z[t][traj.ib1] = traj.b[t]
	end
	return nothing
end

function update_z!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,nq,nu,nw,nc,nb,nz,nθ}
	for t in eachindex((1:traj.H))
		update_z!(traj, t)
	end
	return nothing
end

function update_θ!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}, t::Int) where {T,nq,nu,nw,nc,nb,nz,nθ}
	if t in (1:traj.H)
		traj.θ[t][traj.iq0] = traj.q[t]
		traj.θ[t][traj.iq1] = traj.q[t+1]
		traj.θ[t][traj.iu1] = traj.u[t]
		traj.θ[t][traj.iw1] = traj.w[t]
	end
	return nothing
end

function update_θ!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,nq,nu,nw,nc,nb,nz,nθ}
	for t in eachindex((1:traj.H))
		update_θ!(traj, t)
	end
	return nothing
end
