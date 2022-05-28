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
		traj.u[t] = policy(p, traj, t)

        # disturbances
        traj.w[t] .= disturbances(w, traj.q[t+1], t)

        # step
        status = RoboDojo.step!(s, t, verbose=verbose)
        !status && break
    end

    return status
end
