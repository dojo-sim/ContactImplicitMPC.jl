# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("flamingo", "flat_2D_lc", "flat")
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt));

H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 1000 
κ_mpc = 2.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(model.nu-6); 2; 2]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc])

p = ci_mpc_policy(ref_traj, s, obj,
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
    mpc_opts = CIMPCOptions())

# ## Initial conditions
q1_sim = ContactImplicitMPC.SVector{model.nq}(copy(ref_traj.q[2]))
q0_sim = ContactImplicitMPC.SVector{model.nq}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample))
v1_sim = (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / ref_traj.h

# ## Simulator
sim = Simulator(s, H_sim, h=h_sim, policy=p)

# ## Simulate
simulate!(sim, Array(q1_sim), Array(v1_sim))

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_meshrobot!(vis, model, sim.traj, h=h_sim * 10, sample=10);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim) # Time budget
H_sim * h_sim / sum(sim.stats.dt) # Speed ratio
