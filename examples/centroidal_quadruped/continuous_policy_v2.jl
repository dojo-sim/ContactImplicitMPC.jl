
function simulate!(s::Simulator{T};
	clock_time_noise=1.0e-3,
	verbose=false) where T
    status = false

    p = s.policy
	w = s.dist
	traj = s.traj

	N = length(traj.u)

	# state
	t = 1
	x = [traj.q[t+1]; (traj.q[t+1] - traj.q[t]) ./ s.h; traj.γ[t]]

	# last_computed_control = p.ref_traj.u[1] / N_sample
	clock_time = 0.0
    for t = 1:N
		# println("timestep $t / $N, time $time")
		policy_time = @elapsed traj.u[t] = exec_policy(p, x, clock_time)
		traj.u[t] .*= s.h # transform force from policy to impulse for RoboDojo simulator
		s.opts.record && (s.stats.policy_time[t] = policy_time)

        # disturbances
        traj.w[t] .= disturbances(w, traj.q[t+1], t)

        # step
        status = RoboDojo.step!(s, t, verbose=verbose)
		x .= [traj.q[t+2]; (traj.q[t+2] - traj.q[t+1]) ./ s.h; traj.γ[t]]
		clock_time += s.h + clock_time_noise * rand(1)[1]

        !status && break
    end

    return status
end

function exec_policy(p::CIMPC{T,NQ,NU,NW,NC}, x::Vector{T}, t::T) where {T,NQ,NU,NW,NC}
	# reset
	if t ≈ 0.0
		p.next_time_update = 0.0
		p.q0 .= p.ref_traj.q[1]
		p.altitude .= 0.0
		p.buffer_time = 0.0
		set_trajectory!(p.traj, p.ref_traj)
		set_implicit_trajectory!(p.im_traj, p.im_traj_cache)
		reset_window!(p.window)
	end

    if t >= p.next_time_update
		(p.opts.altitude_update && t > 0.0) && (update_altitude!(p.altitude, p.ϕ, p.s,
									x, NQ, NC, p.N_sample,
									threshold = p.opts.altitude_impact_threshold,
									verbose = p.opts.altitude_verbose))

		set_altitude!(p.im_traj, p.altitude)
		update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)

		# visualize
		# p.opts.live_plotting && live_plotting(p.s.model, p.traj, traj, p.newton, p.q0, traj.q[t+1], t)

		# shift trajectory
		rot_n_stride!(p.traj, p.traj_cache, p.stride, p.window)

		# update
		update_window!(p.window, p.ref_traj.H)

		# reset update time
		p.next_time_update = (t - t % p.traj.h) + p.traj.h
    end

	policy_time = @elapsed begin
		if p.buffer_time <= 0.0
			# @show t
			# optimize
			q1 = x[1:NQ]
			q0 = x[1:NQ] - x[NQ .+ (1:NQ)] .* p.traj.h
			newton_solve!(p.newton, p.s, q0, q1,
				p.window, p.im_traj, p.traj, warm_start = t > 0.0)

			update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)
		end
	end

	# update buffer time
	(p.buffer_time <= 0.0) && (p.buffer_time = policy_time)
	p.buffer_time -= p.traj.h

	# scale control
	p.u .= p.newton.traj.u[1]
	p.u ./= p.traj.h

	return p.u
end
rot_n_stride!
