@testset "Controller: Implicit Dynamics" begin
	# Quadruped reference trajectory
	s = ContactImplicitMPC.get_simulation("quadruped", "flat_2D_lc", "flat")
	model = s.model
	env = s.env

	ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(model, env,
		joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
		load_type = :split_traj_alt))

	ip_opts = InteriorPointOptions(
				κ_tol = 2.0 * 1.0e-4,
				r_tol = 1.0e-8,
				diff_sol = true,
				solver = :empty_solver)

	im_traj = ContactImplicitMPC.ImplicitTrajectory(ref_traj, s,
		κ = 1.0e-4,
		mode = :configuration,
		opts=ip_opts)

	ContactImplicitMPC.implicit_dynamics!(im_traj, ref_traj)

	@test all([norm(im_traj.dq2[t], Inf) < 1.0e-2 for t = 1:ref_traj.H-1])
end
