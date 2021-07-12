@testset "Simulator: Hopper (2D)" begin

    s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    ref_traj = get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
        load_type = :joint_traj)
    T = ref_traj.H
    h = ref_traj.h

    for t = 1:T
    	r = ContactControl.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
    	@test norm(r, Inf) < 1.0e-5
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
    status = ContactControl.simulate!(sim, verbose = false)
    @test status
    @test norm(ref_traj.q[end] - sim.traj.q[end], Inf) < 1.0e-3

end
