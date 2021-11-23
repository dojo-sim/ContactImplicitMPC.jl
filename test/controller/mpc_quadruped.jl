@testset "MPC: Policy for Quadruped" begin
    # ## Simulation
    s = ContactImplicitMPC.get_simulation("quadruped", "flat_2D_lc", "flat");
    model = s.model
    env = s.env

    # ## Reference Trajectory
    ref_traj = deepcopy(get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
        load_type = :split_traj_alt))
    update_friction_coefficient!(ref_traj, model, env);
    
    H = ref_traj.H
    h = ref_traj.h

    # ## MPC setup
    N_sample = 5
    H_mpc = 10
    h_sim = h / N_sample
    H_sim = 1000
    κ_mpc = 2.0e-4

    obj = TrackingObjective(model, env, H_mpc,
        q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.nq-3)]) for t = 1:H_mpc],
        u = [Diagonal(3e-2 * ones(model.nu)) for t = 1:H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

    p = ci_mpc_policy(ref_traj, s, obj,
        H_mpc = H_mpc,
        N_sample = N_sample,
        κ_mpc = κ_mpc,
        mode = :configuration,
        n_opts = NewtonOptions(
            solver = :lu_solver,
            r_tol = 3e-4,
            max_iter = 5,
            ),
        mpc_opts = CIMPCOptions(),
        ip_opts = InteriorPointOptions(
                        undercut = 5.0,
                        γ_reg = 0.1,
                        κ_tol = κ_mpc,
                        r_tol = 1.0e-8,
                        diff_sol = true,
                        solver = :empty_solver,
                        ),
        )

    # ## Initial conditions
    q1_sim, v1_sim = initial_conditions(ref_traj); 

    # ## Simulator
    sim = simulator(s, H_sim, h=h_sim, policy=p);

    # ## Simulate
    status = simulate!(sim, q1_sim, v1_sim);

    @test status 
    qerr, uerr, γerr, berr = ContactImplicitMPC.tracking_error(ref_traj, sim.traj, N_sample, idx_shift=[1])
    @test qerr < 0.0201 * 1.5 # 0.0201
    @test uerr < 0.0437 * 1.5 # 0.0437
    @test γerr < 0.374 * 1.5 # 0.374
    @test berr < 0.0789 * 1.5 # 0.0789
    qerr > 0.0201 * 1.2 && @warn "mild regression on q tracking: current tracking error = $qerr, nominal tracking error = 0.0201"
    uerr > 0.0437 * 1.2 && @warn "mild regression on u tracking: current tracking error = $uerr, nominal tracking error = 0.0437"
    γerr > 0.3740 * 1.2 && @warn "mild regression on γ tracking: current tracking error = $γerr, nominal tracking error = 0.374"
    berr > 0.0789 * 1.2 && @warn "mild regression on b tracking: current tracking error = $berr, nominal tracking error = 0.0789"
end

# @testset "MPC quadruped: long trajectory with structured newton solver" begin
	# s = get_simulation("quadruped", "flat_2D_lc", "flat")
	# model = s.model
	# env = s.env

	# ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
	#     joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
	#     load_type = :split_traj_alt))

	# # time
	# H = ref_traj.H
	# h = ref_traj.h
	# N_sample = 5
	# H_mpc = 10
	# h_sim = h / N_sample
	# H_sim = 1000

	# # barrier parameter
	# κ_mpc = 1.0e-4

	# obj_mpc = ContactImplicitMPC.quadratic_objective(model, H_mpc,
	#     q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.nq-3)]) for t = 1:H_mpc+2],
	#     v = [Diagonal(0.0 * ones(model.nq)) for t = 1:H_mpc],
	#     u = [Diagonal(3e-2 * ones(model.nu)) for t = 1:H_mpc-1])
	# obj_mpc = TrackingObjective(model, env, H_mpc,
	#     q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.nq-3)]) for t = 1:H_mpc],
	#     u = [Diagonal(3e-2 * ones(model.nu)) for t = 1:H_mpc],
	#     γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc],
	#     b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc])

	# p = ci_mpc_policy(ref_traj, s, obj_mpc,
	#     H_mpc = H_mpc,
	#     N_sample = N_sample,
	#     κ_mpc = κ_mpc,
	# 	mode = :configuration,
	# 	newton_mode = :structure,
	#     n_opts = NewtonOptions(
	# 		solver = :lu_solver,
	# 		r_tol = 3e-4,
	# 		max_iter = 5,
	# 		),
	#     mpc_opts = CIMPCOptions(),
	# 	ip_opts = InteriorPointOptions(
	# 		max_iter = 100,
	# 		verbose = false,
	# 		r_tol = 1.0e-4,
	# 		κ_tol = 1.0e-4,
	# 		diff_sol = true,
	# 		solver = :empty_solver,
	# 		),
	#     )
    # # initial configurations
    # q1_sim = ref_traj.q[2]
    # v1_sim = (ref_traj.q[2] - ref_traj.q[1]) ./ h

    # sim = ContactImplicitMPC.simulator(s, H_sim, h=h_sim, policy=p) 

    # # simulator
    # @test status = ContactImplicitMPC.simulate!(sim, q1_sim, v1_sim);
    # ref_traj = deepcopy(ref_traj)
	# qerr, uerr, γerr, berr = tracking_error(ref_traj, sim.traj, N_sample, idx_shift = [1])
	# @test qerr < 0.0202 * 1.5 # 0.0202
	# @test uerr < 0.0437 * 1.5 # 0.0437
	# @test γerr < 0.3780 * 1.5 # 0.3780
	# @test berr < 0.0799 * 1.5 # 0.0799
	# qerr > 0.0201 * 1.2 && @warn "mild regression on q tracking: current tracking error = $qerr, nominal tracking error = 0.0202"
	# uerr > 0.0437 * 1.2 && @warn "mild regression on u tracking: current tracking error = $uerr, nominal tracking error = 0.0437"
	# γerr > 0.3790 * 1.2 && @warn "mild regression on γ tracking: current tracking error = $γerr, nominal tracking error = 0.3780"
	# berr > 0.0798 * 1.2 && @warn "mild regression on b tracking: current tracking error = $berr, nominal tracking error = 0.0799"

# end
