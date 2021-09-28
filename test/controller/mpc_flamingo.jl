@testset "Controller: MPC flamingo" begin
    s = get_simulation("flamingo", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    ref_traj_ = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
        load_type = :split_traj_alt, update_friction = false))
    ref_traj = deepcopy(ref_traj_)

    H = ref_traj.H
    h = ref_traj.h
    N_sample = 5
    H_mpc = 15
    h_sim = h / N_sample
    H_sim = 2100

    # barrier parameter
    κ_mpc = 2.0e-4

    obj = TrackingVelocityObjective(model, env, H_mpc,
        v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
        q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
        u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(model.dim.u-6); 2; 2]) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

    p = linearized_mpc_policy(ref_traj, s, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        n_opts = NewtonOptions(
            r_tol = 3e-4,
            max_iter = 5),
        mpc_opts = LinearizedMPCOptions(
            live_plotting = false,
            altitude_update = false,
            altitude_impact_threshold = 0.02,
            altitude_verbose = false,
            )
        )

    q1_ref = copy(ref_traj.q[2])
    q0_ref = copy(ref_traj.q[1])
    q1_sim = SVector{model.dim.q}(q1_ref)
    q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
    @assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

    sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = InteriorPointOptions(
            undercut = Inf,
            γ_reg = 0.0,
            r_tol = 1.0e-8,
            κ_tol = 1.0e-8),
        sim_opts = SimulatorOptions(warmstart = true)
        )

    @time status = simulate!(sim)
    ref_traj = deepcopy(ref_traj_)

    @test norm(sim.traj.q[1] - ref_traj.q[2]*(1-1/N_sample) - ref_traj.q[1]/N_sample) < 1e-8
    @test norm(sim.traj.q[2] - ref_traj.q[2]) < 1e-8
    @test norm(ref_traj.q[1][2:end] - ref_traj.q[end-1][2:end]) < 1e-8
    @test norm(ref_traj.q[2][2:end] - ref_traj.q[end-0][2:end]) < 1e-8

    qerr, uerr, γerr, berr = ContactImplicitMPC.tracking_error(ref_traj, sim.traj, N_sample, idx_shift=[1])
    @test qerr < 0.0154 * 1.5 # 0.0154
    @test uerr < 0.0829 * 1.5 # 0.0829
    @test γerr < 0.444 * 1.5 # 0.444
    @test berr < 0.0169 * 1.5 # 0.0169
    qerr > 0.0154 * 1.2 && @warn "mild regression on q tracking: current tracking error = $qerr, nominal tracking error = 0.0154"
    uerr > 0.0829 * 1.2 && @warn "mild regression on u tracking: current tracking error = $uerr, nominal tracking error = 0.0829"
    γerr > 0.4440 * 1.2 && @warn "mild regression on γ tracking: current tracking error = $γerr, nominal tracking error = 0.444"
    berr > 0.0169 * 1.2 && @warn "mild regression on b tracking: current tracking error = $berr, nominal tracking error = 0.0169"
end
