@testset "Linearized MPC: Policy for Quadruped" begin
    T = Float64

    s = get_simulation("quadruped", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
        load_type = :split_traj_alt))

    ref_traj = deepcopy(ref_traj_)

    # time
    H = ref_traj.H
    h = ref_traj.h
    N_sample = 5
    H_mpc = 10
    h_sim = h / N_sample
    H_sim = 500

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
    ref_traj = deepcopy(ref_traj_)

    qerr, uerr, γerr, berr = tracking_error(ref_traj, sim.traj, N_sample, idx_shift=[1])
    @test qerr < 0.0201 * 1.5 # 0.0201
    @test uerr < 0.0437 * 1.5 # 0.0437
    @test γerr < 0.374 * 1.5 # 0.374
    @test berr < 0.0789 * 1.5 # 0.0789
    qerr > 0.0201 * 1.2 && @warn "mild regression on q tracking: current tracking error = $qerr, nominal tracking error = 0.0201"
    uerr > 0.0437 * 1.2 && @warn "mild regression on u tracking: current tracking error = $uerr, nominal tracking error = 0.0437"
    γerr > 0.3740 * 1.2 && @warn "mild regression on γ tracking: current tracking error = $γerr, nominal tracking error = 0.374"
    berr > 0.0789 * 1.2 && @warn "mild regression on b tracking: current tracking error = $berr, nominal tracking error = 0.0789"

end

@testset "Linearized MPC: Policy for Quadruped on Sinusoidal Terrain" begin
    T = Float64

    s_sim = get_simulation("quadruped", "sine1_2D_lc", "sinusoidal")
    s = get_simulation("quadruped", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    nq = model.dim.q

    ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
        load_type = :split_traj_alt))
    ref_traj = deepcopy(ref_traj_)

    # time
    H = ref_traj.H
    h = ref_traj.h
    N_sample = 5
    H_mpc = 10
    h_sim = h / N_sample
    H_sim = 1500

    # barrier parameter
    κ_mpc = 1.0e-4

    obj = TrackingObjective(model, env, H_mpc,
        q = [Diagonal(1e-2 * [10; 0.02; 0.25; 0.25 * ones(nq-3)]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

    p = linearized_mpc_policy(ref_traj, s, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        n_opts = NewtonOptions(
            r_tol = 3e-4,
            max_iter = 5,
            # verbose = true,
            ),
        mpc_opts = LinearizedMPCOptions(
            # live_plotting=true,
            altitude_update = true,
            altitude_impact_threshold = 0.05,
            altitude_verbose = false,
            )
        )


    q1_ref = copy(ref_traj.q[2])
    q0_ref = copy(ref_traj.q[1])
    q1_sim = SVector{model.dim.q}(q1_ref)
    q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
    @assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

    sim = simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
        p = p,
        ip_opts = InteriorPointOptions(
            r_tol = 1.0e-8,
            κ_init = 1.0e-6,
            κ_tol = 2.0e-6),
        sim_opts = SimulatorOptions(warmstart = true)
        )

    @time status = simulate!(sim)
    ref_traj = deepcopy(ref_traj_)

    qerr, uerr, γerr, berr = tracking_error(ref_traj, sim.traj, N_sample, idx_shift=[1])
    @test qerr < 0.0333 * 1.5 # 0.0333
    @test uerr < 0.0437 * 1.5 # 0.0437
    @test γerr < 0.381 * 1.5 # 0.381
    @test berr < 0.0795 * 1.5 # 0.0795
    qerr > 0.0333 * 1.2 && @warn "mild regression on q tracking: current tracking error = $qerr, nominal tracking error = 0.0333"
    uerr > 0.0437 * 1.2 && @warn "mild regression on u tracking: current tracking error = $uerr, nominal tracking error = 0.0437"
    γerr > 0.3810 * 1.2 && @warn "mild regression on γ tracking: current tracking error = $γerr, nominal tracking error = 0.381"
    berr > 0.0795 * 1.2 && @warn "mild regression on b tracking: current tracking error = $berr, nominal tracking error = 0.0795"

end
