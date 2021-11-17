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

mutable struct CIMPC{T,NQ,NU,NW,NC,NB,NZ,Nθ,R,RZ,Rθ,Nν,W,FC} <: Policy{T}
	u::Vector{T}
	traj::ContactTraj{T,NQ,NU,NW,NC,NB,NZ,Nθ}
	traj_cache::ContactTraj{T,NQ,NU,NW,NC,NB,NZ,Nθ}
	ref_traj::ContactTraj{T,NQ,NU,NW,NC,NB,NZ,Nθ}
	im_traj::ImplicitTraj{T,R,RZ,Rθ}
	H::Int
	stride::Vector{T}
	altitude::Vector{T}
	ϕ::Vector{T}
	κ::Vector{T}
	newton::Newton{T,NQ,NU,NW,NC,NB,NZ,Nθ,Nν}
	newton_mode::Symbol
	s::Simulation{T,W,FC}
	q0::Vector{T}
	N_sample::Int
	cnt::Vector{Int}
	opts::CIMPCOptions{T}
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
	im_traj = ImplicitTraj(traj, s,
		κ = κ_mpc,
		max_time = mpc_opts.ip_max_time,
		opts=ip_opts,
		mode = mode)
	 
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

	CIMPC(zeros(s.model.nu), traj, traj_cache, ref_traj, im_traj, H_mpc, stride, altitude, ϕ, [κ_mpc], newton, newton_mode, s, copy(ref_traj.q[1]),
		N_sample, [N_sample], mpc_opts)
end


function policy(p::CIMPC, traj, t)
	# reset
	if t == 1
		p.cnt[1] = p.N_sample
		p.q0 .= copy(p.ref_traj.q[1])
	end

    if p.cnt[1] == p.N_sample
		(p.opts.altitude_update && t > 1) && (update_altitude!(p.altitude, p.s,
									traj, t, p.N_sample,
									threshold = p.opts.altitude_impact_threshold,
									verbose = p.opts.altitude_verbose))
		# update!(p.im_traj, p.traj, p.s, p.altitude, κ = p.κ) #@@@ keep the altitude update here
		set_altitude!(p.im_traj, p.altitude) #@@@ keep the altitude update here
		newton_solve!(p.newton, p.s, p.im_traj, p.traj,
			warm_start = t > 1, q0 = copy(p.q0), q1 = copy(traj.q[t+1]))
		update!(p.im_traj, p.traj, p.s, p.altitude, κ = p.κ[1]) #@@@ only keep the rotation stuff not the altitude update.
		p.opts.live_plotting && live_plotting(p.s.model, p.traj, traj, p.newton, p.q0, copy(traj.q[t+1]), t)

		rot_n_stride!(p.traj, p.stride)
		p.q0 .= copy(traj.q[t+1])
		p.cnt[1] = 0
    end

    p.cnt[1] += 1

	if p.newton_mode == :direct
    	return p.newton.traj.u[1] / p.N_sample # rescale output
	elseif p.newton_mode == :structure
		return p.newton.u[1] / p.N_sample
	else
		@error "newton mode specified not available"
	end
end
