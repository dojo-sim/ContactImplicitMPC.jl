"""
	ImplicitTrajectory{T}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linearizedimation computed around `lin`.
"""
mutable struct ImplicitTrajectory{T,R,RZ,Rθ,NQ}
	H::Int
	lin::Vector{LinearizedStep{T}}
	d::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}} # dynamics violation
	dq2::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	dγ1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	db1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	δq0::Vector{SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false}} # q0 solution gradient length=H
	δq1::Vector{SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false}} # q1 solution gradient length=H
	δu1::Vector{SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false}} # u1 solution gradient length=H
	ip::Vector{InteriorPoint{T,R,RZ,Rθ}}
	mode::Symbol
	iq2::SVector{NQ,Int}
end

function ImplicitTrajectory(ref_traj::ContactTraj, s::Simulation;
	κ = ref_traj.κ[1],
	max_time = 1e5,
	mode = :configurationforce,
	opts = InteriorPointOptions(
			undercut = 5.0,
			γ_reg = 0.1,
			κ_tol = κ[1],
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver,
			max_time = max_time))

	model = s.model
	env = s.env

	iq2 = index_q2(model, env, quat = false)
	nq2 = length(iq2)

	H = ref_traj.H

	nq = model.nq
	nu = model.nu
	nc = model.nc
	nb = nc * friction_dim(env)
	if mode == :configurationforce
		nd = nq + nc + nb
	elseif mode == :configuration
		nd = nq
	else
		@error "invalid mode"
	end
	nz = num_var(model, env)
	nθ = num_data(model)

	lin = [LinearizedStep(s, ref_traj.z[t], ref_traj.θ[t], κ) for t = 1:H]

	ip =  [interior_point(
			 deepcopy(ref_traj.z[t]),
			 deepcopy(ref_traj.θ[t]),
			 idx = IndicesOptimization(model, env),
			 r! = rlin!,
			 rz! = rzlin!,
			 rθ! = rθ!,
			 r  = RLin(s, lin[t].z, lin[t].θ, lin[t].r, lin[t].rz, lin[t].rθ),
			 rz = RZLin(s, lin[t].rz),
			 rθ = RθLin(s, lin[t].rθ),
			 opts = opts) for t = 1:H]

	# views
	d = [view(ip[t].z, 1:nd) for t = 1:H]
	dq2 = [view(ip[t].z, 1:nq) for t = 1:H]
	if mode == :configurationforce
		dγ1 = [view(ip[t].z, nq .+ (1:nc)) for t = 1:H]
		db1 = [view(ip[t].z, nq + nc .+ (1:nb)) for t = 1:H]
	elseif mode == :configuration
		dγ1 = [view(ip[t].z, 1:0) for t = 1:H]
		db1 = [view(ip[t].z, 1:0) for t = 1:H]
	else
		@error "invalid mode"
	end

	off = 0
	δq0 = [view(ip[t].δz, 1:nd, off .+ (1:nq)) for t = 1:H]; off += nq
	δq1 = [view(ip[t].δz, 1:nd, off .+ (1:nq)) for t = 1:H]; off += nq
	δu1 = [view(ip[t].δz, 1:nd, off .+ (1:nu)) for t = 1:H]; off += nu

	return ImplicitTrajectory{typeof.([ip[1].z[1], ip[1].r, ip[1].rz, ip[1].rθ])...,nq2}(H, lin, d, dq2, dγ1, db1, δq0, δq1, δu1, ip, mode,
		SVector{nq2,Int}(iq2))
end

function update!(im_traj::ImplicitTrajectory{T,R,RZ,Rθ,nq}, ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
	s::Simulation{T,W,FC}, alt::Vector{T}, κ::T, H::Int) where {T,R,RZ,Rθ,nq,nu,nw,nc,nb,nz,nθ,W,FC}

	for t = 1:H
		# central-path parameter
		im_traj.ip[t].κ[1] = κ

		t > H-1 && continue

		# # residual
		# im_traj.lin[t] = im_traj.lin[t+1]
		# im_traj.ip[t].r = im_traj.ip[t+1].r
		# im_traj.ip[t].rz = im_traj.ip[t+1].rz
		# im_traj.ip[t].rθ = im_traj.ip[t+1].rθ

		# # altitude
		# im_traj.ip[t].r.alt = alt
	end

	# update!(im_traj.lin[H], s, ref_traj.z[H], ref_traj.θ[H])

	# z0  = im_traj.lin[H].z
	# θ0  = im_traj.lin[H].θ
	# r0  = im_traj.lin[H].r
	# rz0 = im_traj.lin[H].rz
	# rθ0 = im_traj.lin[H].rθ

	# update!(im_traj.ip[H].r, z0, θ0, r0, rz0, rθ0)
	# update!(im_traj.ip[H].rz, rz0)
	# update!(im_traj.ip[H].rθ, rθ0)

	# # altitude
	# im_traj.ip[H].r.alt = alt
	return nothing
end

function set_implicit_trajectory!(im_traj::ImplicitTrajectory, im_traj_cache::ImplicitTrajectory)
	H = im_traj.H

	for t = 1:H
		im_traj.lin[t] = im_traj_cache.lin[t]
		im_traj.ip[t].r = im_traj_cache.ip[t].r
		im_traj.ip[t].rz = im_traj_cache.ip[t].rz
		im_traj.ip[t].rθ = im_traj_cache.ip[t].rθ
	end

	return nothing
end

function set_altitude!(im_traj::ImplicitTrajectory, alt::Vector)
	H = im_traj.H
	for t = 1:H
		# altitude
		set_altitude!(im_traj.ip[t], alt)
	end
	return nothing
end

function set_altitude!(ip::InteriorPoint, alt::Vector)
	# altitude
	ip.r.alt = alt
	return nothing
end

function implicit_dynamics!(im_traj::ImplicitTrajectory, traj::ContactTraj;
	threads=false,
	window=collect(1:traj.H + 2))

	for (i, t) in enumerate(window[1:end-2])
		# initialized solver
		z_initialize!(im_traj.ip[t].z, im_traj.iq2, traj.q[i+2]) #TODO: try alt. schemes
		im_traj.ip[t].θ .= traj.θ[i]
	end

	if threads
		Threads.@threads for t in window[1:end-2]
			# solve
			status = interior_point_solve!(im_traj.ip[t])
			!status && (@warn "implicit dynamics failure (t = $t)")
		end
	else
		for t in window[1:end-2]
			# solve
			status = interior_point_solve!(im_traj.ip[t])
			!status && (@warn "implicit dynamics failure (t = $t)")
		end
	end

	for (i, t) in enumerate(window[1:end-2])
		# compute dynamics violation
		im_traj.dq2[t] .-= traj.q[i+2]

		if im_traj.mode == :configurationforce
			im_traj.dγ1[t] .-= traj.γ[i]
			im_traj.db1[t] .-= traj.b[i]
		elseif im_traj.mode == :configuration
			nothing
		end
	end
	return nothing
end
