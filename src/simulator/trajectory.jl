struct ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}
	H::Int
	h::T
	κ::Vector{T}
	q::Vector#{SizedArray{Tuple{nq},T,1,1}}         # trajectory of q's   length=H+2
	u::Vector#{SizedArray{Tuple{nu},T,1,1}}         # trajectory of u's   length=H
	w::Vector#{SizedArray{Tuple{nw},T,1,1}}         # trajectory of w's   length=H
	γ::Vector#{SizedArray{Tuple{nc},T,1,1}}         # trajectory of γ's   length=H
	b::Vector#{SizedArray{Tuple{nb},T,1,1}}         # trajectory of b's   length=H
	z::Vector#{SizedArray{Tuple{nz},T,1,1}}         # trajectory of z's   length=H
	θ::Vector#{SizedArray{Tuple{nθ},T,1,1}}         # trajectory of θ's   length=H
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

	# q = [zeros(SizedVector{nq,T}) for k=1:H+2]
	# u = [zeros(SizedVector{nu,T}) for k=1:H]
	# w = [zeros(SizedVector{nw,T}) for k=1:H]
	# γ = [zeros(SizedVector{nc,T}) for k=1:H]
	# b = [zeros(SizedVector{nb,T}) for k=1:H]
	# z = [zeros(SizedVector{nz,T}) for k=1:H]
	# θ = [zeros(SizedVector{nθ,T}) for k=1:H]

	q = [zeros(nq) for k=1:H+2]
	u = [zeros(nu) for k=1:H]
	w = [zeros(nw) for k=1:H]
	γ = [zeros(nc) for k=1:H]
	b = [zeros(nb) for k=1:H]
	z = [zeros(nz) for k=1:H]
	θ = [[zeros(nθ-1); copy(h)] for k=1:H]
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

function repeat_ref_traj(traj::ContactTraj, model::ContactDynamicsModel, N::Int;
       idx_shift = (1:0))

    shift = (traj.q[end] - traj.q[2])[idx_shift]

    q = [deepcopy(traj.q)...]
    u = [deepcopy(traj.u)...]
    w = [deepcopy(traj.w)...]
    γ = [deepcopy(traj.γ)...]
    b = [deepcopy(traj.b)...]
    z = [deepcopy(traj.z)...]
    θ = [deepcopy(traj.θ)...]
    κ = deepcopy(traj.κ)

    for i = 1:N-1
		for t = 1:traj.H
			qt = copy(traj.q[t + 2])
			qt[idx_shift] .+= i * shift
			push!(q, qt)
			push!(u, copy(traj.u[t]))
			push!(w, copy(traj.w[t]))
			push!(γ, copy(traj.γ[t]))
			push!(b, copy(traj.b[t]))
			push!(z, copy(traj.z[t]))
			push!(θ, copy(traj.θ[t]))
		end
    end

    return ContactTraj{typeof(traj).parameters...}(traj.H * N, traj.h, traj.κ,
        q, u, w, γ, b, z, θ,
		traj.iq0, traj.iq1, traj.iu1, traj.iw1, traj.iq2, traj.iγ1, traj.ib1)
end

function sub_ref_traj(traj::ContactTraj, model::ContactDynamicsModel, idx)

    q = deepcopy(traj.q[collect([[idx..., idx[1][end] + 1, idx[1][end] + 2]]...)])
    u = deepcopy(traj.u[idx])
    w = deepcopy(traj.w[idx])
    γ = deepcopy(traj.γ[idx])
    b = deepcopy(traj.b[idx])
    z = deepcopy(traj.z[idx])
    θ = deepcopy(traj.θ[idx])
    κ = deepcopy(traj.κ)

    return ContactTraj{typeof(traj).parameters...}(length(idx), traj.h, traj.κ,
        q, u, w, γ, b, z, θ,
		traj.iq0, traj.iq1, traj.iu1, traj.iw1, traj.iq2, traj.iγ1, traj.ib1)
end

function update_friction_coefficient!(traj::ContactTraj, model::ContactDynamicsModel)
	for t = 1:traj.H
		q2, γ1, b1, ψ1, η1, __ = unpack_z(model, traj.z[t])
		q0, q1, u1, w1, _, h = unpack_θ(model, traj.θ[t])
		traj.z[t] .= pack_z(model, q2, γ1, b1, ψ1, η1)
		traj.θ[t] .= pack_θ(model, q0, q1, u1, w1, copy(model.μ_world), h)
	end
	nothing
end
