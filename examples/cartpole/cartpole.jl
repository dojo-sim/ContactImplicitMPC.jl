# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("walledcartpole", "flat_2D_lc", "flat");
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = contact_trajectory(model, env, H, h);
qref = [0.0; 0.0; 0.0; 0.0]
ur = zeros(model.nu)
γr = zeros(model.nc)
br = zeros(model.nc * friction_dim(env))
ψr = zeros(model.nc)
ηr = zeros(model.nc * friction_dim(env))
wr = zeros(model.nw)
h = 0.04
H = 50

# ## Set Reference
for t = 1:H
	ref_traj.z[t] = pack_z(model, env, qref, γr, br, ψr, ηr)
	ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
end

# ## Simulation Time
h_sim = h / N_sample
H_sim = 1000

# ## No Policy
p = open_loop_policy(fill(zeros(nu), H_sim); N_sample = N_sample)

# ## MPC setup 
N_sample = 2
H_mpc = 10
κ_mpc = 2.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [[Diagonal([1e-1;1e-3;1e-8;1e-8]) for t = 1:H_mpc-1];
		 [Diagonal([1e+1;1e+0;1e-8;1e-8]) for t = 1:1]],
	v = [[Diagonal([1e-0;3e+1;1e-8;1e-8]) for t = 1:H_mpc-1];
		 [Diagonal([1e+1;1e+1;1e-8;1e-8]) for t = 1:1]],
	u = [Diagonal([3e-2]) for t = 1:H_mpc])

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    mpc_opts =CIMPCOptions(),
	n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		));

# ## Disturbances
idx_d1 = 20
idx_d2 = idx_d1 + 200
idx_d3 = idx_d2 + 150
idx_d4 = idx_d3 + 200
idx_d5 = idx_d4 + 150
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
impulses = [[+0.2; 0; 0; 0], [+0.2; 0; 0; 0], [-0.2; 0; 0; 0], [+0.2; 0; 0; 0], [+0.2; 0; 0; 0]]
d = impulse_disturbances(impulses, idx);

# ## Initial Conditions
q1_sim = zeros(4) 
v1_sim = zeros(4) 

# ## Simulator
sim = simulator(s, H_sim, h=h_sim, policy=p, dist=d)

# ## Simulate
status = simulate!(sim, q1_sim, v1_sim);;

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_robot!(vis, model, sim.traj, h=h_sim, sample = 1)

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio