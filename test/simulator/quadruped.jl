@testset "Simulator: Quadruped" begin
    # Reference trajectory
	s = get_simulation("quadruped", "flat_2D_lc", "flat")
	model = s.model
	env = s.env
    model.μ_world = 0.5

    ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
		joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
		load_type = :split_traj_alt))
    ContactControl.update_friction_coefficient!(ref_traj, model, env)

    T = ref_traj.H
    h = ref_traj.h

    for t = 1:T
    	r = ContactControl.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
    	@test norm(r) < 1.0e-4
    end

    # initial conditions
    q0 = SVector{model.dim.q}(ref_traj.q[1])
    q1 = SVector{model.dim.q}(ref_traj.q[2])

    # simulator
    sim = ContactControl.simulator(s, q0, q1, h, T,
        p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u]),
        ip_opts = ContactControl.InteriorPointOptions(
    		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-5, solver = :lu_solver),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    # simulate
    @test status = ContactControl.simulate!(sim, verbose = false)
    # @show sim.traj.q[end][1:3]
    @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.025
end
