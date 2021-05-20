@testset "Controller: MPC flamingo" begin
    s = get_simulation("flamingo", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
        joinpath(pwd(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
        load_type = :split_traj_alt))

    H = ref_traj.H
    h = ref_traj.h
    N_sample = 5
    H_mpc = 15
    h_sim = h / N_sample
    H_sim = 1000#35000

    # barrier parameter
    κ_mpc = 1.0e-4

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
            # live_plotting=true,
            # altitude_update = true,
            altitude_impact_threshold = 0.02,
            altitude_verbose = true,
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
            r_tol = 1.0e-8,
            κ_init = 1.0e-8,
            κ_tol = 2.0e-8),
        sim_opts = SimulatorOptions(warmstart = true)
        )

    @time status = simulate!(sim)
end
