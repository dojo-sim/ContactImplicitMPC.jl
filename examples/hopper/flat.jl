# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays

# ## Simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj));

H = ref_traj.H
h = ref_traj.h

# ## MPC setup 
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 1000
κ_mpc = 2.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [0.1,3,1,3])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-3, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc]);

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
	ip_opts = InteriorPointOptions(
					undercut = 5.0,
					κ_tol = κ_mpc,
					r_tol = 1.0e-8,
					diff_sol = true,
					solver = :empty_solver,
					max_time = 1e5),
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions());

# ## Initial conditions
q1_sim = SVector{model.dim.q}(copy(ref_traj.q[2]))
q0_sim = SVector{model.dim.q}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample))

# ## Simulator
sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
		γ_reg = 0.0,
		undercut = Inf,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8,),
    sim_opts = SimulatorOptions(warmstart = true));

# ## Simulate
status = simulate!(sim, verbose = true);

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
open(vis)

# ## Visualize
anim = visualize_robot!(vis, model, sim.traj, sample=5);
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample = 5);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim) # Time budget
H_sim * h_sim / sum(sim.stats.dt) # Speed ratio



