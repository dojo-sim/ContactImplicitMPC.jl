# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("quadruped", "piecewise1_2D_lc", "piecewise", approx = true);
model = s.model
env = s.env
nq = model.dim.q
nc = model.dim.c

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt));

H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 4000 #3000
κ_mpc = 1.0e-4

obj = TrackingObjective(model, env, H_mpc,
	q = [Diagonal(1e-2 * [5; 0.02; 0.10; 0.25 * ones(nq-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc]);

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
		solver = :lu_solver,
		r_tol = 3e-4,
		max_iter = 5,
		),
    mpc_opts = CIMPCOptions(
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        ),
	ip_opts = InteriorPointOptions(
		max_iter = 100,
		verbose = false,
		r_tol = 1.0e-4,
		κ_tol = 1.0e-4,
		diff_sol = true,
		solver = :empty_solver,
		),
    )

# ## Initial conditions
q1_sim = ContactImplicitMPC.SVector{model.dim.q}(copy(ref_traj.q[2]))
q0_sim = ContactImplicitMPC.SVector{model.dim.q}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample));

# ## Simulator
sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	ip_opts = InteriorPointOptions(
		undercut = Inf,
		γ_reg = 0.0,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = SimulatorOptions(warmstart = true),
    );

# ## Simulation
@time status = simulate!(sim, verbose = true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
plot_surface!(vis, s.env, ylims=[0.3, -0.05])
anim = visualize_meshrobot!(vis, model, sim.traj, sample=5);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim) # Time budget
H_sim * h_sim / sum(sim.stats.dt) # Speed ratio