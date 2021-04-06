"""
    linearized model-predictive control policy
"""
mutable struct LinearizedMPC <: Policy
   ref_traj
   im_traj
   stride
   altitude
   κ
   newton
   model
   q0
   N_sample
   cnt
   live_plotting
end

function linearized_mpc_policy(traj, model, cost;
	H_mpc = ref_traj.H,
	N_sample = 1,
	κ_mpc = traj.κ[1],
	n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		verbose = false,
		live_plotting = false),
	ip_max_time = 100.0,
	live_plotting = false)

	ref_traj = deepcopy(traj)

	im_traj = ImplicitTraj(ref_traj, model, κ = κ_mpc, max_time = ip_max_time)
	stride = get_stride(model, ref_traj)
	altitude = zeros(model.dim.c)
	newton = Newton(H_mpc, ref_traj.h, model, cost = cost, opts = n_opts)

	LinearizedMPC(ref_traj, im_traj, stride, altitude, κ_mpc, newton, model, copy(ref_traj.q[1]),
		N_sample, N_sample, live_plotting)
end


function policy(p::LinearizedMPC, x, traj, t)
    if p.cnt == p.N_sample
		update_altitude!(p.altitude, x)
		update!(p.im_traj, p.ref_traj, p.model, p.altitude, κ = p.κ)
		newton_solve!(p.newton, p.model, p.im_traj, p.ref_traj,
			verbose = p.newton.opts.verbose, warm_start = t > 0.0, q0 = copy(p.q0), q1 = copy(x))
		rot_n_stride!(p.ref_traj, p.stride)
		p.q0 .= copy(x)
		p.cnt = 0

		p.live_plotting && live_plotting(p.model, p.traj, traj, p.newton)
    end

    p.cnt += 1

    return p.newton.traj.u[1] / p.N_sample # rescale output
end
