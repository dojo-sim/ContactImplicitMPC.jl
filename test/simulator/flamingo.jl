@testset "Simulator: Flamingo" begin
    # Reference trajectory
	s = get_simulation("flamingo", "flat_2D_lc", "flat")
	model = s.model
	env = s.env

	model.μ_world = 0.1
	ref_traj = get_trajectory(s.model, s.env,
		joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait1.jld2"),
		load_type = :split_traj_alt)
	ContactImplicitMPC.update_friction_coefficient!(ref_traj, model, env)

	for t = 1:ref_traj.H
		r = ContactImplicitMPC.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
		@test norm(r) < 1.0e-4
	end

	# initial conditions
	q0 = SVector{model.nq}(ref_traj.q[1])
	q1 = SVector{model.nq}(ref_traj.q[2])

	h = ref_traj.h
	H = ref_traj.H
	model.μ_world = 0.9

	"""
	    PD tracking policy
	"""
	mutable struct PDPolicy <: ContactImplicitMPC.Policy
		model::Model
	    traj::ContactTrajectory
		q0::AbstractVector
	    idx::Int
	    cnt::Int
	    N_sample::Int
	end

	function pd_policy(model, traj; N_sample = 1)
	    PDPolicy(model, traj, copy(traj.q[1]), 0, N_sample, N_sample)
	end

	function ContactImplicitMPC.policy(p::PDPolicy, x, traj, t)
	    # reset
	    if t == 1
	        p.idx = 0
	        p.cnt = p.N_sample
	    end

	    if p.cnt == p.N_sample
	        p.idx += 1
	        p.cnt = 0
			p.q0 .= copy(x)
	    end

	    p.cnt += 1

		u = p.traj.u[p.idx]
		# PD
		kp = 1.0 * ones(p.model.nu)
		kp[1] *= 100.0
		# kd = 10.0

		# u = Diagonal(kp) * B_func(p.model, x) * (x - p.traj.q[p.idx + 1])
		# @show u
		# u -= kd * B_func(p.model, x) * (x - traj.q[t]) / traj.h

	    return u ./ p.N_sample
	end

	p = pd_policy(model, ref_traj)

	# simulator
	sim = ContactImplicitMPC.simulator(s, q0, q1, h, 24,
	    p = p,
	    ip_opts = ContactImplicitMPC.InteriorPointOptions(
			r_tol = 1.0e-8,
			κ_tol = 1.0e-8),
	    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true))

	# simulate
	status = ContactImplicitMPC.simulate!(sim, verbose = false)
	@test status
end

# include(joinpath(module_dir(), "src", "dynamics", "flamingo", "visuals.jl"))
# vis = Visualizer()
# open(vis)
# anim = visualize_robot!(vis, model, sim.traj, name = :sim)
# anim = visualize_robot!(vis, model, ref_traj, anim = anim, name = :ref)
