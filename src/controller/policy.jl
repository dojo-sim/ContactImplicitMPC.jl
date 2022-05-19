"""
    contact-implicit model-predictive control policy
"""

@with_kw mutable struct CIMPCOptions{T}
	altitude_update::Bool = false
	altitude_impact_threshold::T = 1.0
	altitude_verbose::Bool = false
    ip_max_time::T = 1e5     # maximum time allowed for an InteriorPoint solve
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct CIMPC{T,NQ,NU,NW,NC,NB,NZ,Nθ,R,RZ,Rθ,Nν,W,FC,NQQ,NJ,NR,NI,OB,LS,NV} <: Policy{T}
	u::Vector{T}
	traj::ContactTraj{T,NQ,NU,NW,NC,NB,NZ,Nθ}
	traj_cache::ContactTraj{T,NQ,NU,NW,NC,NB,NZ,Nθ}
	ref_traj::ContactTraj{T,NQ,NU,NW,NC,NB,NZ,Nθ}
	im_traj::ImplicitTrajectory{T,R,RZ,Rθ,NQQ}
	im_traj_cache::ImplicitTrajectory{T,R,RZ,Rθ,NQQ}
	H::Int
	stride::Vector{T}
	altitude::Vector{T}
	ϕ::Vector{T}
	κ::Vector{T}
	newton::Newton{T,NQ,NU,NW,NC,NB,NZ,Nθ,Nν,NJ,NR,NI,OB,LS,NV}
	newton_mode::Symbol
	s::Simulation{T,W,FC}
	q0::Vector{T}
	N_sample::Int
	cnt::Vector{Int}
	window::Vector{Int}
	opts::CIMPCOptions{T}

	buffer_time::T
	next_time_update::T
end

function ci_mpc_policy(traj::ContactTraj, s::Simulation{T}, obj::Objective;
	H_mpc = traj.H,
	N_sample = 1,
	κ_mpc = traj.κ[1],
	mode = :configurationforce,
	newton_mode = :direct,
	n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		verbose = false,
		live_plotting = false),
	mpc_opts = CIMPCOptions(),
	ip_opts = InteriorPointOptions(
				γ_reg = 0.1,
				undercut = 5.0,
				κ_tol = κ_mpc,
				r_tol = 1.0e-8,
				diff_sol = true,
				solver = :empty_solver,
				max_time = mpc_opts.ip_max_time,)
	) where T

	traj = deepcopy(traj)
	traj_cache = deepcopy(traj)
	ref_traj = deepcopy(traj)
	im_traj = ImplicitTrajectory(traj, s,
		κ = κ_mpc,
		max_time = mpc_opts.ip_max_time,
		opts=ip_opts,
		mode = mode)

	im_traj_cache = deepcopy(im_traj)

	stride = get_stride(s.model, traj)
	altitude = zeros(s.model.nc)
	ϕ = zeros(s.model.nc)
	if newton_mode == :direct
		newton = Newton(s, H_mpc, traj.h, traj, im_traj, obj = obj, opts = n_opts)
	elseif newton_mode == :structure
		newton = NewtonStructure(s, H_mpc, traj, obj, κ_mpc, opts = n_opts)
	else
		@error "invalid Newton solver specified"
	end

	window = zeros(Int, H_mpc + 2)


	CIMPC(zeros(s.model.nu), traj, traj_cache, ref_traj, im_traj, im_traj_cache,
		H_mpc, stride, altitude, ϕ, [κ_mpc], newton, newton_mode, s, copy(ref_traj.q[1]),
		N_sample, [N_sample], window, mpc_opts,
		0.0, 0.0)
end

function policy(p::CIMPC{T,NQ,NU,NW,NC}, traj::Trajectory{T}, t::Int) where {T,NQ,NU,NW,NC}
	# reset
	if t == 1
		p.cnt[1] = p.N_sample
		p.q0 .= p.ref_traj.q[1]
		p.altitude .= 0.0
		set_trajectory!(p.traj, p.ref_traj)
		set_implicit_trajectory!(p.im_traj, p.im_traj_cache)
		reset_window!(p.window)
	end

    if p.cnt[1] == p.N_sample
		# altitude_update
		(p.opts.altitude_update && t > 1) && (update_altitude!(p.altitude, p.ϕ, p.s,
									traj, t, NC, p.N_sample,
									threshold = p.opts.altitude_impact_threshold,
									verbose = p.opts.altitude_verbose))
		set_altitude!(p.im_traj, p.altitude)

		# # optimize
		q1 = traj.q[t+1]
		newton_solve!(p.newton, p.s, p.q0, q1,
			p.window, p.im_traj, p.traj, warm_start = t > 1)

		update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)

		# visualize
		p.opts.live_plotting && live_plotting(p.s.model, p.traj, traj, p.newton, p.q0, traj.q[t+1], t)

		# shift trajectory
		rot_n_stride!(p.traj, p.traj_cache, p.stride, p.window)

		# update
		update_window!(p.window, p.ref_traj.H)
		p.q0 .= q1

		# reset count
		p.cnt[1] = 0
    end

    p.cnt[1] += 1

	# scale control
	if p.newton_mode == :direct
		p.u .= p.newton.traj.u[1]
		p.u ./= p.N_sample
	elseif p.newton_mode == :structure
		p.u .= p.newton.u[1]
		p.u ./= p.N_sample
	else
		println("newton mode specified not available")
	end

	return p.u
end

function reset_window!(window)
	n = length(window)
	for i = 1:n
		window[i] = i
	end
	return
end

function update_window!(window, max_window)
	n = length(window)
	for i = 1:n
		window[i] += 1
		if window[i] > max_window
			window[i] = 1
		end
	end
	return window
end
