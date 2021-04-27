@testset "Linearized MPC: Policy for Quadruped" begin
    T = Float64

    # get model
    model = ContactControl.get_model("quadruped")

    # get trajectory
    ref_traj = get_trajectory("quadruped", "gait2", load_type = :split_traj_alt)
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
    obj = TrackingObjective(H_mpc, model.dim,
        q = [Diagonal(1e-2 * [10; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

    # linearized MPC policy
    p = ContactControl.linearized_mpc_policy(ref_traj, model, obj,
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
    sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_init = 1.0e-8,
            κ_tol = 2.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    # simulator
    @test status = ContactControl.simulate!(sim)

    # # test lengths
    # @test p.newton.traj.H == H_mpc
    # @test length(p.newton.traj.q) == H_mpc+2
    # @test length(p.newton.traj.u) == H_mpc
    # @test p.newton.traj_cand.H == H_mpc
    # @test length(p.newton.traj.q) == H_mpc+2
    # @test length(p.newton.traj.u) == H_mpc
    # @test length(p.newton.ν) == H_mpc
    # @test length(p.newton.ν_cand) == H_mpc
    #
    # # Check tracking performance
    # function tracking_error(ref_traj::ContactControl.ContactTraj, sim::ContactControl.Simulator, p::ContactControl.Policy, M)
    #     N_sample = p.N_sample
    #     # M = p.idx
    #     qq = []
    #     for q in ref_traj_copy.q
    #         for i = 1:N_sample
    #             push!(qq, q)
    #         end
    #     end
    #     uu = []
    #     γγ = []
    #     bb = []
    #     for t = 1:ref_traj_copy.H
    #         for i = 1:N_sample
    #             push!(uu, ref_traj_copy.u[t])
    #             push!(γγ, ref_traj_copy.γ[t])
    #             push!(bb, ref_traj_copy.b[t])
    #         end
    #     end
    #
    #     q_error = abs(mean(hcat(qq...)[1:model.dim.q, 1:M] - hcat(sim.traj.q...)[1:model.dim.q, 1:M]))
    #     u_error = abs(mean(hcat(uu...)[1:model.dim.u, 1:M] - hcat(sim.traj.u...)[1:model.dim.u, 1:M]))
    #     γ_error = abs(mean(hcat(γγ...)[1:model.dim.c, 1:M] - hcat(sim.traj.γ...)[1:model.dim.c, 1:M]))
    #     b_error = abs(mean(hcat(bb...)[1:model.dim.b, 1:M] - hcat(sim.traj.b...)[1:model.dim.b, 1:M]))
    #
    #     return q_error, u_error, γ_error, b_error
    # end
    #
    # # Check the tracking error with disturbances
    # q_error, u_error, γ_error, b_error = tracking_error(ref_traj, sim, p, 100)
    #
    # @test q_error < 0.10
    # @test u_error < 0.10
    # @test γ_error < 0.19
    # @test b_error < 0.15
end

@testset "Linearized MPC: Policy for Quadruped on Sinusoidal Terrain" begin
    T = Float64
    # get hopper model
    model_sim = get_model("quadruped", surf="sinusoidal")

    model = get_model("quadruped", surf="flat")
    nq = model.dim.q
    nu = model.dim.u
    nc = model.dim.c
    nb = model.dim.b
    nd = nq + nc + nb
    nr = nq + nu + nc + nb + nd

    # get trajectory
    ref_traj = get_trajectory("quadruped", "gait2", load_type = :split_traj_alt)
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

    obj = TrackingObjective(H_mpc, model.dim,
        q = [Diagonal(1e-2 * [0.02; 0.02; 1.0; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

    p = linearized_mpc_policy(ref_traj, model, obj,
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
    q1_sim = SVector{model.dim.q}(q1_ref)
    q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
    @assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

    w_amp = [+0.02, -0.20]
    sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        # d = open_loop_disturbances([rand(model.dim.w) .* w_amp for i=1:H_sim]),
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

        q_error = abs(mean(hcat(qq...)[1:model.dim.q, 1:M] - hcat(sim.traj.q...)[1:model.dim.q, 1:M]))
        u_error = abs(mean(hcat(uu...)[1:model.dim.u, 1:M] - hcat(sim.traj.u...)[1:model.dim.u, 1:M]))
        γ_error = abs(mean(hcat(γγ...)[1:model.dim.c, 1:M] - hcat(sim.traj.γ...)[1:model.dim.c, 1:M]))
        b_error = abs(mean(hcat(bb...)[1:model.dim.b, 1:M] - hcat(sim.traj.b...)[1:model.dim.b, 1:M]))

        return q_error, u_error, γ_error, b_error
    end

    # Check the tracking error with disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, sim, p, 100)

    @test q_error < 0.10
    @test u_error < 0.10
    @test γ_error < 0.19
    @test b_error < 0.15
end

@testset "Linearized MPC: Policy for Particle on Sinusoidal Path" begin
    T = Float64

    # get model
    model = get_model("particle", surf = "sinusoidal")

    # get ref. trajectory
    ref_traj = ContactControl.get_trajectory("particle", "sinusoidal2",
        model_name = "particle_sinusoidal")

    ref_traj_copy = deepcopy(ref_traj)

    # time
    H = ref_traj.H
    h = ref_traj.h

    # initial conditions
    q0 = SVector{model.dim.q}(ref_traj.q[1])
    q1 = SVector{model.dim.q}(ref_traj.q[2])

    # simulator
    sim = ContactControl.simulator(model, q0, q1, 1.0 * h, H,
        p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u], N_sample = 1),
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_init = 1.0e-8, κ_tol = 2.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    @time status = ContactControl.simulate!(sim)

    # linearized motion planning
    obj = ContactControl.TrackingObjective(H, model.dim,
        q = [Diagonal(1.0 * ones(model.dim.q))    for t = 1:H],
        u = [Diagonal(1.0e-1 * ones(model.dim.u)) for t = 1:H],
        γ = [Diagonal(1.0e-6 * ones(model.dim.c)) for t = 1:H],
        b = [Diagonal(1.0e-6 * ones(model.dim.b)) for t = 1:H])

    # Simulate linearized MPC policy
    κ_mpc = 1.0e-4
    tf = h * H
    H_mpc = 20
    N_sample = 5
    h_sim = h / N_sample
    H_sim = 500
    q1_ref = copy(ref_traj.q[2])
    q0_ref = copy(copy(ref_traj.q[1]))
    q1_sim = SVector{model.dim.q}(q1_ref)
    q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
    @assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

    p = ContactControl.linearized_mpc_policy(ref_traj, model, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        n_opts = ContactControl.NewtonOptions(
            r_tol = 3e-4,
            max_iter = 5),
        mpc_opts = ContactControl.LinearizedMPCOptions())

    # simulator
    sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_init = 1.0e-8,
            κ_tol = 2.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    @time status = ContactControl.simulate!(sim)

    # test lengths
    @test p.newton.traj.H == H_mpc
    @test length(p.newton.traj.q) == H_mpc+2
    @test length(p.newton.traj.u) == H_mpc
    @test p.newton.traj_cand.H == H_mpc
    @test length(p.newton.traj.q) == H_mpc+2
    @test length(p.newton.traj.u) == H_mpc
    @test length(p.newton.ν) == H_mpc
    @test length(p.newton.ν_cand) == H_mpc

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

        q_error = abs(mean(hcat(qq...)[1:model.dim.q, 1:M] - hcat(sim.traj.q...)[1:model.dim.q, 1:M]))
        u_error = abs(mean(hcat(uu...)[1:model.dim.u, 1:M] - hcat(sim.traj.u...)[1:model.dim.u, 1:M]))
        γ_error = abs(mean(hcat(γγ...)[1:model.dim.c, 1:M] - hcat(sim.traj.γ...)[1:model.dim.c, 1:M]))
        b_error = abs(mean(hcat(bb...)[1:model.dim.b, 1:M] - hcat(sim.traj.b...)[1:model.dim.b, 1:M]))

        return q_error, u_error, γ_error, b_error
    end

    # Check the tracking error with disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, sim, p, 100)

    @test q_error < 0.01
    @test u_error < 0.01
    @test γ_error < 0.001
    @test b_error < 0.001
end
