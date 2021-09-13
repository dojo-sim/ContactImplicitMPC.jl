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

function contact_trajectory(model::ContactModel, env::Environment, H::Int, h::T; κ::T=0.0) where {T}
	dim = model.dim
	nq = dim.q
    nu = dim.u
    nw = dim.w
    nc = dim.c
    nb = nc * friction_dim(env)
	nz = num_var(model, env)
	nθ = num_data(model)

	q = [zeros(nq) for k=1:H+2]
	u = [zeros(nu) for k=1:H]
	w = [zeros(nw) for k=1:H]
	γ = [zeros(nc) for k=1:H]
	b = [zeros(nb) for k=1:H]
	z = [zeros(nz) for k=1:H]
	θ = [[zeros(nθ-1); copy(h)] for k=1:H]
	κ = [κ]
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
	#TODO: specify type
	vqq
	vqqq
	vqu
	vγq
	vγqq
	vγu
	vbq
	vbqq
	vbu
end

function contact_derivative_trajectory(model::ContactModel, env::Environment, δz::AbstractArray, H::Int)
	nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = nc * friction_dim(env)

	dq2dq0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:H]
	dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:H]
	dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:H]
	dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:H]
	dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:H]
	dγdu = [SizedMatrix{nc,nu}(zeros(nc, nu)) for t = 1:H]
	dbdq0 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:H]
	dbdq1 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:H]
	dbdu = [SizedMatrix{nb,nu}(zeros(nb, nu)) for t = 1:H]

	vqq = view(δz, 1:nq, 1:nq)
	vqqq = view(δz, 1:nq, nq .+ (1:nq))
	vqu = view(δz, 1:nq, 2 * nq .+ (1:nu))
	vγq = view(δz, nq .+ (1:nc), 1:nq)
	vγqq = view(δz, nq .+ (1:nc), nq .+ (1:nq))
	vγu = view(δz, nq .+ (1:nc), 2 * nq .+ (1:nu))
	vbq = view(δz, nq + nc .+ (1:nb), 1:nq)
	vbqq = view(δz, nq + nc .+ (1:nb), nq .+ (1:nq))
	vbu = view(δz, nq + nc .+ (1:nb), 2 * nq .+ (1:nu))

	ContactDerivTraj(
		dq2dq0,
		dq2dq1,
		dq2du,
		dγdq0,
		dγdq1,
		dγdu,
		dbdq0,
		dbdq1,
		dbdu,
		vqq,
		vqqq,
		vqu,
		vγq,
		vγqq,
		vγu,
		vbq,
		vbqq,
		vbu)
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

function repeat_ref_traj(traj::ContactTraj, N::Int;
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

function sub_traj(traj::ContactTraj, idx::AbstractVector{Int})

    q = deepcopy(traj.q[collect([[idx..., idx[end] + 1, idx[end] + 2]]...)])
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

function update_friction_coefficient!(traj::ContactTraj, model::ContactModel, env::Environment)
	for t = 1:traj.H
		q2, γ1, b1, ψ1, s1, η1, __ = unpack_z(model, env, traj.z[t])
		q0, q1, u1, w1, _, h = unpack_θ(model, traj.θ[t])
		traj.z[t] .= pack_z(model, env, q2, γ1, b1, ψ1, η1)
		traj.θ[t] .= pack_θ(model, q0, q1, u1, w1, copy(model.μ_world), h)
	end
	nothing
end

function get_trajectory(model::ContactModel, env::Environment, gait_path::String;
		load_type::Symbol = :split_traj, update_friction::Bool = false)
	#TODO: assert model exists

	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = nc * friction_dim(env)
	res = JLD2.jldopen(gait_path)

	if load_type == :split_traj
		q, u, γ, b, h = res["q"], res["u"], res["γ"], res["b"], mean(res["h̄"])
		ū = res["ū"]
		ψ = [ut[nu + nc + nb .+ (1:nc)] for ut in ū]
		η = [ut[nu + nc + nb + nc .+ (1:nb)] for ut in ū]

		T = length(u)

		traj = contact_trajectory(T, h, model, env)
		traj.q .= deepcopy(q)
		traj.u .= deepcopy(u)
		traj.γ .= deepcopy(γ)
		traj.b .= deepcopy(b)
		traj.z .= [pack_z(model, env, q[t+2], γ[t], b[t], ψ[t], η[t]) for t = 1:T]
		traj.θ .= [pack_θ(model, q[t], q[t+1], u[t], zeros(nw), model.μ_world, h) for t = 1:T]
	elseif load_type == :split_traj_alt
		q, u, γ, b, ψ, η, μ, h = res["qm"], res["um"], res["γm"], res["bm"], res["ψm"], res["ηm"], res["μm"], res["hm"]
		T = length(u)

		traj = contact_trajectory(model, env, T, h)
		traj.q .= deepcopy(q)
		traj.u .= deepcopy(u)
		traj.γ .= deepcopy(γ)
		traj.b .= deepcopy(b)
		traj.z .= [pack_z(model, env, q[t+2], γ[t], b[t], ψ[t], η[t]) for t = 1:T]
		traj.θ .= [pack_θ(model, q[t], q[t+1], u[t], zeros(nw), μ, h) for t = 1:T]
	elseif load_type == :joint_traj
		traj = res["traj"]
	end
	update_friction && update_friction_coefficient!(traj, model, env)
	return traj
end

# Check tracking performance
function tracking_error(ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
		sim_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}, N_sample::Int; idx_shift=(1:0)) where {T,nq,nu,nw,nc,nb,nz,nθ}

	# Horizons
	H_ref = ref_traj.H
	H_sim = sim_traj.H
	H̄_sim = H_sim / N_sample
	dupl_traj = repeat_ref_traj(ref_traj, Int(ceil(H̄_sim / H_ref)), idx_shift=idx_shift)
	H_dupl = dupl_traj.H

	q_error = 0.0
	u_error = 0.0
	γ_error = 0.0
	b_error = 0.0
	cnt = 0
    for t = 1:H_dupl
		cnt += 1
		(t-1)*N_sample + 1 > H_sim && break
		q_error += norm(dupl_traj.q[t+2] - sim_traj.q[(t-1)*N_sample + 3], 1) / nq
		u_error += norm(dupl_traj.u[t] - sim_traj.u[(t-1)*N_sample + 1], 1) / nu
		γ_error += norm(dupl_traj.γ[t] - sim_traj.γ[(t-1)*N_sample + 1], 1) / nc
		b_error += norm(dupl_traj.b[t] - sim_traj.b[(t-1)*N_sample + 1], 1) / nb
    end

	q_error /= cnt
	u_error /= cnt
	γ_error /= cnt
	b_error /= cnt
    return q_error, u_error, γ_error, b_error
end
