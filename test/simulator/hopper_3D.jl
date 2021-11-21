@testset "Simulator: Hopper (3D)" begin
    # Reference trajectory
    s = get_simulation("hopper_3D", "flat_3D_lc", "flat")
    model = s.model
    env = s.env

    ref_traj = get_trajectory(s.model, s.env,
    	joinpath(module_dir(), "src/dynamics/hopper_3D/gaits/gait_in_place.jld2"),
    	load_type = :joint_traj)
    T = ref_traj.H
    h = ref_traj.h

    for t = 1:T
    	r = ContactImplicitMPC.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
    	@test norm(r, Inf) < 1.0e-5
    end

    # initial conditions
    q1 = ref_traj.q[2]
    v1 = (ref_traj.q[2] - ref_traj.q[1]) ./ h

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h)

    # simulate
    status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test status
end
