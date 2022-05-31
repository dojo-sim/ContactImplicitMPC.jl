
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
		# p.next_time_update = 0.0
		p.q0 .= p.ref_traj.q[1]
		# p.altitude .= 0.0
		p.buffer_time = 0.0
		# set_trajectory!(p.traj, p.ref_traj)
		set_implicit_trajectory!(p.im_traj, p.im_traj_cache)
		reset_window!(p.window)
	end

	policy_time = @elapsed begin
		if p.buffer_time <= 0.0

			# window 
			idx_nearest = findnearest(p.times_reference, t % (p.ref_traj.h * (p.ref_traj.H - 1)))[1]
			for i = 1:(p.H + 2)
				p.window[i] = i + idx_nearest - 1 > p.ref_traj.H ? i - p.ref_traj.H + idx_nearest - 1 : i + idx_nearest - 1
			end
			set_window!(p.traj, p.ref_traj, p.window)

			# optimize
			q1 = x[1:NQ]
			q0 = x[1:NQ] - x[NQ .+ (1:NQ)] .* p.traj.h
			newton_solve!(p.newton, p.s, q0, q1,
				p.window, p.im_traj, p.traj, warm_start = t > 0.0)
		end
	end

	# update buffer time
	(p.buffer_time <= 0.0) && (p.buffer_time = policy_time)
	p.buffer_time -= p.traj.h

	# control
	p.u .= p.newton.traj.u[1]

	# add gains
	if p.opts.gains
		K1 = p.K_traj[p.window[1]]

		q1 = x[1:NQ]
		q0 = x[1:NQ] - x[NQ .+ (1:NQ)] .* p.traj.h
		p.u .+= K1 * ([p.traj.q[p.window[1]]; p.traj.q[p.window[2]]] - [q0; q1])
	end
	
	p.u ./= p.traj.h

	return p.u
end


# horizon_reference = 11
# horizon_mpc = 3
# timestep = 0.1
# time_reference = timestep * (horizon_reference - 1)
# times_reference = range(0, stop=time_reference, length=horizon_reference)

# t = time_reference - 0.1

function findnearest(a,x)
	length(a) > 0 || return 0:-1
	r = searchsorted(a,x)
	length(r) > 0 && return r
	last(r) < 1 && return searchsorted(a,a[first(r)])
	first(r) > length(a) && return searchsorted(a,a[last(r)])
	x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
	x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
	return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end

# idx_nearest = findnearest(times_reference, t % time_reference)[1]
# idx_window_wrap = [i + idx_nearest - 1 > horizon_reference ? i - horizon_reference + idx_nearest - 1 : i + idx_nearest - 1 for i in 1:(horizon_mpc + 2)]



# ###
# # window 
# t = 0.2
# p.window .= 0
# idx_nearest = findnearest(p.times_reference, t % (p.ref_traj.h * (p.ref_traj.H - 1)))[1]
# for i = 1:(p.H + 2)
# 	p.window[i] = i + idx_nearest - 1 > p.ref_traj.H ? i - p.ref_traj.H + idx_nearest - 1 : i + idx_nearest - 1
# end
# p.window
# @show idx_nearest
# @show length(p.traj.q)
# @show length(p.traj.u) 
# @show length(p.ref_traj.q) 
# @show length(p.ref_traj.u) 
# @show length(p.window) 
# @show p.window
# set_window!(p.traj, p.ref_traj, p.window)