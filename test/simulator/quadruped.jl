@testset "Simulator: Quadruped" begin
    # Reference trajectory
    s = get_simulation("quadruped", "flat_2D_lc", "flat")
    model = s.model
    env = s.env
    model.μ_world = 0.5

    ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
		joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
		load_type = :split_traj_alt))
    ContactImplicitMPC.update_friction_coefficient!(ref_traj, model, env)

    T = ref_traj.H
    h = ref_traj.h

    for t = 1:T
      r = ContactImplicitMPC.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
      @test norm(r) < 1.0e-4
    end

    # initial conditions
    v1 = (ref_traj.q[2] - ref_traj.q[1]) ./ h
    q1 = ref_traj.q[2]

    # policy 
    p = ContactImplicitMPC.open_loop_policy([ut for ut in ref_traj.u]) 

    # simulator
    sim = ContactImplicitMPC.simulator(s, T, h=h, policy=p)

    # simulate
    @test status = ContactImplicitMPC.simulate!(sim, q1, v1)
    @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.025
end
