# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays

# ## Simulation
s_sim = get_simulation("hopper_3D", "flat_3D_lc", "flat");
model_sim = s_sim.model
env_sim = s_sim.env

s = get_simulation("hopper_3D", "flat_3D_lc", "flat")
model = s.model
env = s.env

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_3D/gaits/gait_forward.jld2"),
    load_type = :joint_traj));

H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 10
H_mpc = 20
h_sim = h / N_sample
H_sim = 1200 # 12000
κ_mpc = 1.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [3,3,0.1,5e+1,5e+1,5e+1,10])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-1, 1e-1, 1e1]) for t = 1:H_mpc],
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
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        )
    )

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		max_time = ref_traj.h, # HARD REAL TIME
		),
    mpc_opts = LinearizedMPCOptions(
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        ),
	ip_opts = InteriorPointOptions(
		max_iter = 100,
		verbose = false,
		r_tol = 1.0e-4,
		κ_tol = 1.0e-4,
		diff_sol = true,
		solver = :empty_solver,
		),
    );

# ## Initial conditions
q1_sim = SVector{model.dim.q}(copy(ref_traj.q[2]))
q0_sim = SVector{model.dim.q}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample))

# ## Simulator
sim = simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        γ_reg = 0.0,
        undercut = Inf,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = SimulatorOptions(warmstart = true));

# ## Simulate
@time status = ContactImplicitMPC.simulate!(sim, verbose = true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
open(vis)

# ## Visualize
plot_surface!(vis, s_sim.env, n=200)
visualize_robot!(vis, model, sim.traj)

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim) # Time budget
H_sim * h_sim / sum(sim.stats.dt) # Speed ratio
