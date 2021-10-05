# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s_no_load = get_simulation("quadruped", "flat_2D_lc", "flat");
s_load = get_simulation("quadruped", "flat_2D_lc", "payload",
    model_variable_name="quadruped_payload",
    dynamics_name="dynamics_payload");

model_no_load = s_no_load.model
model_load = s_load.model
env_no_load = s_no_load.env
env_load = s_load.env

# ## Reference Trajectory
ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s_no_load.model, s_no_load.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt));
H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 2100
κ_mpc = 1.0e-4

obj = TrackingObjective(model_no_load, env_no_load, H_mpc,
    q = [Diagonal(1e-2 * [10.0; 0.02; 0.25; 0.5 * ones(model_no_load.dim.q-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model_no_load.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model_no_load.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model_no_load.dim.c * friction_dim(env_no_load))) for t = 1:H_mpc]);

p = ci_mpc_policy(ref_traj, s_no_load, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        verbose = false,
        max_iter = 5),
    mpc_opts = CIMPCOptions(
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        )
    );

# ## Initial conditions
q1_sim = ContactImplicitMPC.SVector{model_load.dim.q}(copy(ref_traj.q[2]))
q0_sim = ContactImplicitMPC.SVector{model_load.dim.q}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample))

# ## Simulator
sim_load = simulator(s_load, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        undercut = Inf,
		γ_reg = 0.0,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = SimulatorOptions(warmstart=true));

# ## Simulate
@time status = ContactImplicitMPC.simulate!(sim_load, verbose=true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_meshrobot!(vis, s_load.model, sim_load.traj, sample=1, name=:Payload);
anim = visualize_payload!(vis, s_load.model, sim_load.traj, anim=anim, sample=1, name=:Payload, object=:mesh);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim) # Time budget
H_sim * h_sim / sum(sim.stats.dt) # Speed ratio
