
function simulate!(s::Simulator{T}; verbose=false) where T
    status = false

    p = s.policy
	w = s.dist
	traj = s.traj
	N = length(traj.u)

	# last_computed_control = p.ref_traj.u[1] / N_sample
	buffer_time = 0.0
    for t = 1:N
		# println("timestep $t / $N, time $time")
		policy_time = @elapsed traj.u[t] = exec_policy(p, traj, t, buffer_time)
		s.opts.record && (s.stats.policy_time[t] = policy_time)
		(buffer_time < 0.0) && (buffer_time = policy_time)
		buffer_time -= s.h # decrement the buffer time of the sim step

		# traj.u[t] .= last_computed_control
		# if buffer_time <= 0.0
			# policy
			# buffer_time = @elapsed last_computed_control = exec_policy(p, traj, t, buffer_time)
			# s.opts.record && (s.stats.policy_time[t] = buffer_time)
        # end
		# buffer_time -= s.h # decrement the buffer time of the sim step

        # disturbances
        traj.w[t] .= disturbances(w, traj.q[t+1], t)

        # step
        status = RoboDojo.step!(s, t, verbose=verbose)
        !status && break
    end

    return status
end


function exec_policy(p::CIMPC{T,NQ,NU,NW,NC}, traj::Trajectory{T}, t::Int, buffer_time::T) where {T,NQ,NU,NW,NC}
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

		# # # optimize
		# q1 = traj.q[t+1]
		# newton_solve!(p.newton, p.s, p.q0, q1,
		# 	p.window, p.im_traj, p.traj, warm_start = t > 1)

		update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)

		# visualize
		p.opts.live_plotting && live_plotting(p.s.model, p.traj, traj, p.newton, p.q0, traj.q[t+1], t)

		# shift trajectory
		rot_n_stride!(p.traj, p.traj_cache, p.stride, p.window)

		# update
		update_window!(p.window, p.ref_traj.H)
		# p.q0 .= q1

		# reset count
		p.cnt[1] = 0
    end
	p.cnt[1] += 1


	if buffer_time <= 0.0
		@show t
		# optimize
		q0 = traj.q[max(1,t+1-N_sample)]
		q1 = traj.q[t+1]
		newton_solve!(p.newton, p.s, q0, q1,
			p.window, p.im_traj, p.traj, warm_start = t > 1)

		update!(p.im_traj, p.traj, p.s, p.altitude, p.κ[1], p.traj.H)
	end


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
