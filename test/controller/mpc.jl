@testset "Linearized MPC: dummy MPC" begin
    T = Float64
    # get model
    model = ContactControl.get_model("quadruped")

    # get trajectory

    ref_traj = ContactControl.get_trajectory("quadruped", "gait1", load_type=:split_traj)
    H = ref_traj.H
    h = ref_traj.h
    κ = 1.0e-4

    n_opts = ContactControl.NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
    m_opts = ContactControl.MPCOptions{T}(
                N_sample=10,
                M=2*H,
                H_mpc=10,
                κ=κ,
                κ_sim=1e-8,
                r_tol_sim=1e-8,
                open_loop_mpc=false,
                w_amp=[0.05, 0.00],
                live_plotting=false)

    cost = ContactControl.CostFunction(H, model.dim,
        q = [Diagonal(1e-2 * [0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,]) for t = 1:m_opts.H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:m_opts.H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:m_opts.H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:m_opts.H_mpc])
    core = ContactControl.Newton(m_opts.H_mpc, h, model, cost=cost, opts=n_opts)
    mpc = ContactControl.MPC(model, ref_traj, m_opts=m_opts)
    ContactControl.dummy_mpc(model, core, mpc)

    # Check length
    N_sample = mpc.m_opts.N_sample
    M = mpc.m_opts.M
    H_mpc = mpc.m_opts.H_mpc

    @test length(mpc.q_sim) == N_sample*M+2
    @test length(mpc.u_sim) == N_sample*M
    @test length(mpc.w_sim) == N_sample*M
    @test length(mpc.γ_sim) == N_sample*M
    @test length(mpc.b_sim) == N_sample*M

    @test core.traj.H == H_mpc
    @test length(core.traj.q) == H_mpc+2
    @test length(core.traj.u) == H_mpc
    @test core.traj_cand.H == H_mpc
    @test length(core.traj.q) == H_mpc+2
    @test length(core.traj.u) == H_mpc
    @test length(core.ν) == H_mpc
    @test length(core.ν_cand) == H_mpc

    # Check tracking performance

    function tracking_error(ref_traj::ContactControl.ContactTraj, mpc::ContactControl.MPC)
        N_sample = mpc.m_opts.N_sample
        M = mpc.m_opts.M
        q_error = []
        u_error = []
        γ_error = []
        b_error = []
        q_sim = mpc.q_sim[3:N_sample:end]
        u_sim = mpc.u_sim[1:N_sample:end]
        γ_sim = mpc.γ_sim[1:N_sample:end]
        b_sim = mpc.b_sim[1:N_sample:end]
        for t = 1:M
            push!(q_error, norm(ref_traj.q[3+(t-1)%H][2:end] - q_sim[t][2:end]))
            push!(u_error, norm(ref_traj.u[1+(t-1)%H] - u_sim[t]*N_sample))
            push!(γ_error, norm(ref_traj.γ[1+(t-1)%H] - γ_sim[t]*N_sample))
            push!(b_error, norm(ref_traj.b[1+(t-1)%H] - b_sim[t]*N_sample))

        end
        return mean(q_error), mean(u_error), mean(γ_error), mean(b_error)
    end

    # Check the tracking error with disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, mpc)

    @test q_error < 0.10
    @test u_error < 0.10
    @test γ_error < 0.19
    @test b_error < 0.15


    # get trajectory
    ref_traj = ContactControl.get_trajectory("quadruped", "gait1", load_type=:split_traj)
    H = ref_traj.H
    h = ref_traj.h
    κ = 1.0e-4

    n_opts = ContactControl.NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
    m_opts = ContactControl.MPCOptions{T}(
                N_sample=10,
                M=2*H,
                H_mpc=10,
                κ=κ,
                κ_sim=1e-8,
                r_tol_sim=1e-8,
                open_loop_mpc=false,
                w_amp=zeros(model.dim.w),
                live_plotting=false)
    cost = ContactControl.CostFunction(H, model.dim,
        q = [Diagonal(1e-2 * [0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,]) for t = 1:m_opts.H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:m_opts.H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:m_opts.H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:m_opts.H_mpc])
    core = ContactControl.Newton(m_opts.H_mpc, h, model, cost=cost, opts=n_opts)
    mpc = ContactControl.MPC(model, ref_traj, m_opts=m_opts)
    ContactControl.dummy_mpc(model, core, mpc)

    # Check the tracking error with no disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, mpc)

    @test q_error < 0.10
    @test u_error < 0.01
    @test γ_error < 0.15
    @test b_error < 0.10
end

@testset "Linearized MPC: Policy for Quadruped" begin
    T = Float64

    # get model
    model = ContactControl.get_model("quadruped")

    # get trajectory
    ref_traj = ContactControl.get_trajectory("quadruped", "gait1")
    ref_traj_copy = deepcopy(ref_traj)

    # time
    H = ref_traj.H
    h = ref_traj.h
    N_sample = 2
    H_mpc = 10
    h_sim = h / N_sample
    H_sim = 200

    # barrier parameter
    κ = 1.0e-4

    # cost
    cost = ContactControl.CostFunction(H_mpc, model.dim,
        q = [Diagonal(1e-2 * [0.02, 0.02, 1.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

    # linearized MPC policy
    p = ContactControl.linearized_mpc_policy(ref_traj, model, cost, H_mpc, 0, N_sample,
        n_opts = ContactControl.NewtonOptions(r_tol = 3e-4,
            κ_init = κ,
            κ_tol = 2κ,
            solver_inner_iter = 5),
        m_opts = ContactControl.MPCOptions{Float64}(
                    N_sample = N_sample,
                    M = 200,
                    H_mpc = H_mpc,
                    κ = κ,
                    κ_sim = 1e-8,
                    r_tol_sim = 1e-8,
                    open_loop_mpc = false,
                    w_amp = [-0.10, -0.10],
                    ip_max_time = 0.1,
                    live_plotting = false))

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
    @time status = ContactControl.simulate!(sim)


    # test lengths
    @test p.core.traj.H == H_mpc
    @test length(p.core.traj.q) == H_mpc+2
    @test length(p.core.traj.u) == H_mpc
    @test p.core.traj_cand.H == H_mpc
    @test length(p.core.traj.q) == H_mpc+2
    @test length(p.core.traj.u) == H_mpc
    @test length(p.core.ν) == H_mpc
    @test length(p.core.ν_cand) == H_mpc

    # Check tracking performance
    function tracking_error(ref_traj::ContactControl.ContactTraj, sim::ContactControl.Simulator, p::ContactControl.Policy, M)
        N_sample = p.mpc.m_opts.N_sample
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

    ref_traj.κ[1] = 1.0e-4
    ref_traj_copy = deepcopy(ref_traj)

    # time
    H = ref_traj.H
    h = ref_traj.h

    # initial conditions
    q0 = SVector{model.dim.q}(ref_traj.q[1])
    q1 = SVector{model.dim.q}(ref_traj.q[2])

    # simulator
    sim = ContactControl.simulator(model, q0, q1, 1.0 * h, H,
        p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u], h, N_sample = 1),
        ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_init = 1.0e-8, κ_tol = 2.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = false))

    # simulate
    @time status = ContactControl.simulate!(sim)

    # plot(hcat(ref_traj.q...)[1:3, :]',
    #     label = ["x" "y" "z"], color = :black, width = 3.0)
    # plot!(hcat(sim.traj.q...)[1:3, :]',
    #     label = ["x" "y" "z"], color = :red, width = 1.0, legend = :topleft)

    # linearized motion planning
    cost = ContactControl.CostFunction(H, model.dim,
        q = [Diagonal(1.0 * ones(model.dim.q))    for t = 1:H],
        u = [Diagonal(1.0e-1 * ones(model.dim.u)) for t = 1:H],
        γ = [Diagonal(1.0e-6 * ones(model.dim.c)) for t = 1:H],
        b = [Diagonal(1.0e-6 * ones(model.dim.b)) for t = 1:H])

    # Simulate linearized MPC policy
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

    p = ContactControl.linearized_mpc_policy(ref_traj, model, cost, H_mpc, 0, N_sample)

    # # simulator
    sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = ContactControl.InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_init = 1.0e-8,
            κ_tol = 2.0e-8),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    @time status = ContactControl.simulate!(sim)

    # test lengths
    @test p.core.traj.H == H_mpc
    @test length(p.core.traj.q) == H_mpc+2
    @test length(p.core.traj.u) == H_mpc
    @test p.core.traj_cand.H == H_mpc
    @test length(p.core.traj.q) == H_mpc+2
    @test length(p.core.traj.u) == H_mpc
    @test length(p.core.ν) == H_mpc
    @test length(p.core.ν_cand) == H_mpc

    # Check tracking performance
    function tracking_error(ref_traj::ContactControl.ContactTraj, sim::ContactControl.Simulator, p::ContactControl.Policy, M)
        N_sample = p.mpc.m_opts.N_sample
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
