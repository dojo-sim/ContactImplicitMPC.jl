# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("quadruped", "flat_2D_lc", "flat");
model = s.model
env = s.env

model.μ_world = 0.35
# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/calipso_gait11.jld2"),
    load_type = :split_traj_alt))
update_friction_coefficient!(ref_traj, model, env);
H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 1500
κ_mpc = 2.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.nq-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-3 * ones(model.nu)) for t = 1:H_mpc],
	v = [Diagonal(1.0e-5 * ones(model.nq)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
    n_opts = NewtonOptions(
		solver = :ldl_solver,
		r_tol = 3e-4,
		max_iter = 20,
		),
    mpc_opts = CIMPCOptions(),
	ip_opts = InteriorPointOptions(
					undercut = 2.0,
					γ_reg = 0.1,
					κ_tol = κ_mpc,
					r_tol = 1.0e-6,
					diff_sol = true,
					solver = :empty_solver,
					# max_time = 1000.0,
					),
    )

# ## Initial conditions
q1_sim, v1_sim = initial_conditions(ref_traj); 

# ## Simulator
sim = simulator(s, H_sim, h=h_sim, policy=p);

# ## Simulate
simulate!(sim, q1_sim, v1_sim)

# ## Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.open(vis)

# ## Visualize
anim = visualize_meshrobot!(vis, model, sim.traj, h=h_sim * 5, sample=5);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio