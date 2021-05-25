@testset "Linearized MPC: Policy for Quadruped" begin
    T = Float64

    s = get_simulation("quadruped", "flat_2D_lc", "flat")
   	model = s.model
   	env = s.env
   	# s.model.μ_world = 0.5

   	ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
   		joinpath(pwd(), "src/dynamics/quadruped/gaits/gait2.jld2"),
   		load_type = :split_traj_alt))
   	# ContactControl.update_friction_coefficient!(ref_traj, s.model, s.env)

    ref_traj_copy = deepcopy(ref_traj)

    # time
    H = ref_traj.H
    h = ref_traj.h
    N_sample = 2
    H_mpc = 10
    h_sim = h / N_sample
    H_sim = 200

    # barrier parameter
    κ_mpc = 1.0e-4

    # obj
    obj = TrackingVelocityObjective(model, env, H_mpc,
        q = [Diagonal(1e-2 * [10; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
        v = [Diagonal(0e-2 * [10; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

    # linearized MPC policy
    p = ContactControl.linearized_mpc_policy(ref_traj, s, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        n_opts = ContactControl.NewtonOptions(
            r_tol = 3e-4,
            max_iter = 5),
        mpc_opts = ContactControl.LinearizedMPCOptions())

    # initial configurations
    q1_ref = copy(ref_traj.q[2])
    q0_ref = copy(copy(ref_traj.q[1]))
    q1_sim = SVector{model.dim.q}(q1_ref)
    q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
    @assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

    # simulator
    sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_init = 1.0e-8,
            κ_tol = 2.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    # simulator
    @test status = ContactControl.simulate!(sim)
end

@testset "Linearized MPC: Policy for Quadruped on Sinusoidal Terrain" begin
    T = Float64
    # get hopper model
    s_sim = get_simulation("quadruped", "sine1_2D_lc", "flat")
    s_mpc = get_simulation("quadruped", "flat_2D_lc", "flat")

    model = s_sim.model
    env = s_sim.env

    nq = model.dim.q
    nu = model.dim.u
    nc = model.dim.c
    nb = nc * friction_dim(env)
    nd = nq + nc + nb
    nr = nq + nu + nc + nb + nd

    # get trajectory
    ref_traj = deepcopy(ContactControl.get_trajectory(s_mpc.model, s_mpc.env,
        joinpath(pwd(), "src/dynamics/quadruped/gaits/gait2.jld2"),
        load_type = :split_traj_alt))
    ref_traj_copy = deepcopy(ref_traj)

    # time
    H = ref_traj.H
    h = ref_traj.h
    N_sample = 5
    H_mpc = 10
    h_sim = h / N_sample
    H_sim = 500

    # barrier parameter
    κ_mpc = 1.0e-4

    obj = TrackingObjective(s_mpc.model, s_mpc.env, H_mpc,
        q = [Diagonal(1e-2 * [0.02; 0.02; 1.0; 0.25 * ones(s_mpc.model.dim.q-3)]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(s_mpc.model.dim.u)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(s_mpc.model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(s_mpc.model.dim.c * friction_dim(s_mpc.env))) for t = 1:H_mpc])

    p = linearized_mpc_policy(ref_traj, s_mpc, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        n_opts = NewtonOptions(
            r_tol = 3e-4,
            max_iter = 5),
        mpc_opts = LinearizedMPCOptions(
            # live_plotting=true,
            altitude_update = false,
            altitude_impact_threshold = 0.05,
            altitude_verbose = false,
            )
        )

    q1_ref = copy(ref_traj.q[2])
    q0_ref = copy(ref_traj.q[1])
    q1_sim = SVector{s_mpc.model.dim.q}(q1_ref)
    q0_sim = SVector{s_mpc.model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
    @test norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

    sim = simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_init = 1.0e-6,
            κ_tol = 2.0e-6),
        sim_opts = SimulatorOptions(warmstart = true)
        )

    @test status = simulate!(sim)

    # Check tracking performance
    function tracking_error(ref_traj::ContactControl.ContactTraj, sim::ContactControl.Simulator, p::ContactControl.Policy, M)
        N_sample = p.N_sample
        # M = p.idx
        qq = []
        for q in ref_traj_copy.q
            for i = 1:N_sample
                push!(qq, q)
            end
        end
        uu = []
        γγ = []
        bb = []
        for t = 1:ref_traj_copy.H
            for i = 1:N_sample
                push!(uu, ref_traj_copy.u[t])
                push!(γγ, ref_traj_copy.γ[t])
                push!(bb, ref_traj_copy.b[t])
            end
        end

        q_error = abs(mean(hcat(qq...)[1:s_mpc.model.dim.q, 1:M] - hcat(sim.traj.q...)[1:s_mpc.model.dim.q, 1:M]))
        u_error = abs(mean(hcat(uu...)[1:s_mpc.model.dim.u, 1:M] - hcat(sim.traj.u...)[1:s_mpc.model.dim.u, 1:M]))
        γ_error = abs(mean(hcat(γγ...)[1:s_mpc.model.dim.c, 1:M] - hcat(sim.traj.γ...)[1:s_mpc.model.dim.c, 1:M]))
        b_error = abs(mean(hcat(bb...)[1:s_mpc.model.dim.c * friction_dim(s_mpc.env), 1:M] - hcat(sim.traj.b...)[1:s_mpc.model.dim.c * friction_dim(s_mpc.env), 1:M]))

        return q_error, u_error, γ_error, b_error
    end

    # Check the tracking error with disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, sim, p, 100)

    @test q_error < 0.10
    @test u_error < 0.10
    # @test γ_error < 0.19
    # @test b_error < 0.15
end
