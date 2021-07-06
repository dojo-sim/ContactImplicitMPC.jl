"""
	ImplicitTraj{T}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linearizedimation computed around `lin`.
"""
mutable struct ImplicitTraj{T}
	H::Int
	lin::Vector{LinearizedStep{T}}
	d::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}} # dynamics violation
	dq2::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	dγ1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	db1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	# δq0::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # q0 solution gradient length=H
	# δq1::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # q1 solution gradient length=H
	# δu1::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # u1 solution gradient length=H
	δq0::Vector{SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}} # q0 solution gradient length=H
	δq1::Vector{SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}}  # q1 solution gradient length=H
	δu1::Vector{SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}}  # u1 solution gradient length=H
	ip::Vector{<:AbstractIPSolver}
	mode::Symbol
end


function ImplicitTraj(ref_traj::ContactTraj, s::Simulation;
	κ = ref_traj.κ[1],
	max_time = 60.0,
	mode = :configurationforce,
	ip_type::Symbol = :interior_point,
	opts = eval(interior_point_options(ip_type))(
			κ_init = κ[1],
			κ_tol = 2.0 * κ[1],
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver,
			max_time = max_time))

	model = s.model
	env = s.env

	H = ref_traj.H

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
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


	@warn "different init of ips"
	ip =  [eval(ip_type)(
			 # zeros(num_var(model, env)),
			 # zeros(num_data(model)),
			 deepcopy(ref_traj.z[t]),
			 deepcopy(ref_traj.θ[t]),
			 idx_ineq = inequality_indices(model, env),
			 ix = linearization_var_index(model, env)[1],
			 iy1 = linearization_var_index(model, env)[2],
			 iy2 = linearization_var_index(model, env)[3],
			 ibil = linearization_term_index(model, env)[3],
			 r! = r!,
			 rm! = rm!,
			 rz! = rz!,
			 rθ! = rθ!,
			 r  = RLin(s, lin[t].z, lin[t].θ, lin[t].r, lin[t].rz, lin[t].rθ),
			 # rm  = RLin(s, lin[t].z, lin[t].θ, lin[t].r, lin[t].rz, lin[t].rθ),
			 rz = RZLin(s, lin[t].rz),
			 rθ = RθLin(s, lin[t].rθ),
			 v_pr = view(zeros(1,1), 1,1),
			 v_du = view(zeros(1,1), 1,1),
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

	return ImplicitTraj(H, lin, d, dq2, dγ1, db1, δq0, δq1, δu1, ip, mode)
end


function update!(im_traj::ImplicitTraj, ref_traj::ContactTraj,
	s::Simulation,
	 alt::Vector
	; κ = ref_traj.κ[1])

	H = ref_traj.H
	for t = 1:H
		fill!(im_traj.ip[t].κ, κ)
	end
	for t = 1:H-1
		im_traj.lin[t]   = im_traj.lin[t+1]
		im_traj.ip[t].r  = im_traj.ip[t+1].r
		# @warn "need to remove comment"
		im_traj.ip[t].r̄  = im_traj.ip[t+1].r̄
		im_traj.ip[t].rz = im_traj.ip[t+1].rz
		im_traj.ip[t].rθ = im_traj.ip[t+1].rθ

		# altitude
		im_traj.ip[t].r.alt = alt
		# @warn "need to remove comment"
		im_traj.ip[t].r̄.alt = alt
	end

	update!(im_traj.lin[H], s, ref_traj.z[H], ref_traj.θ[H])

	z0  = im_traj.lin[H].z
	θ0  = im_traj.lin[H].θ
	r0  = im_traj.lin[H].r
	rz0 = im_traj.lin[H].rz
	rθ0 = im_traj.lin[H].rθ

	update!(im_traj.ip[H].r,  z0, θ0, r0, rz0, rθ0)
	# @warn "need to remove comment"
	update!(im_traj.ip[H].r̄,  z0, θ0, r0, rz0, rθ0)
	update!(im_traj.ip[H].rz, rz0)
	update!(im_traj.ip[H].rθ, rθ0)

	# altitude
	im_traj.ip[H].r.alt = alt
	# @warn "need to remove comment"
	im_traj.ip[H].r̄.alt = alt
	return nothing
end

function set_altitude!(im_traj::ImplicitTraj, alt::Vector)
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
	ip.r̄.alt = alt
	return nothing
end

function set_altitude!(ip::Mehrotra, alt::Vector)
	# altitude
	ip.r.alt = alt
	return nothing
end

"""
	implicit_dynamics!(im_traj::ImplicitTraj, model::ContactModel,
		traj::ContactTraj; κ = traj.κ)
Compute the evaluations and Jacobians of the implicit dynamics on the trajectory 'traj'. The computation is
linearized since it relies on a linearization about a reference trajectory.
"""
function implicit_dynamics!(im_traj::ImplicitTraj, s::Simulation,
	traj; κ = traj.κ) where {T, D}

	model = s.model
	env = s.env

	# Threads.@threads for t = 1:traj.H
	for t = 1:traj.H
		@show traj.H
		# initialized solver
		z_initialize!(im_traj.ip[t].z, model, env, copy(traj.q[t+2])) #TODO: try alt. schemes
		im_traj.ip[t].θ .= traj.θ[t]

		# solve
		status = interior_point_solve!(im_traj.ip[t])

		# !status && error("implicit dynamics failure (t = $t)")
		!status && (@warn "implicit dynamics failure (t = $t)")

		# compute dynamics violation
		# println("before!: im_traj.dq2[t] = ", scn.(im_traj.dq2[t]))
		im_traj.dq2[t] .-= traj.q[t+2]
		# println("implicit_dynamics!:    traj.q[t+2] = ", scn.(traj.q[t+2]))
		# println(" after!: im_traj.dq2[t] = ", scn.(im_traj.dq2[t]))

		if im_traj.mode == :configurationforce
			im_traj.dγ1[t] .-= traj.γ[t]
			im_traj.db1[t] .-= traj.b[t]
		elseif im_traj.mode == :configuration
			nothing
		end
	end
	return nothing
end
