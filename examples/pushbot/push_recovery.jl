# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("pushbot", "flat_2D_lc", "flat");
model = s.model
env = s.env

# ## Reference Trajectory
h = 0.04
H = 100
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h
qref = [0.0; 0.0]
ur = zeros(model.nu)
γr = zeros(model.nc)
br = zeros(model.nc * friction_dim(env))
ψr = zeros(model.nc)
ηr = zeros(model.nc * friction_dim(env))
wr = zeros(model.nw)

# ## Set Reference
for t = 1:H
	ref_traj.z[t] = pack_z(model, env, qref, γr, br, ψr, ηr)
	ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
end

# ## Initial conditions
# q0 = ContactImplicitMPC.SVector{2}([0.0 * π, 0.0])
q1 = [0.0 * π, 0.0]
v1 = [0.0; 0.0]

# ## Simulator
sim = simulator(s, H, h=h)

# ## Simulate
status = simulate!(sim, q1, v1)

# ## MPC setup 
N_sample = 2
H_mpc = 40
h_sim = h / N_sample
H_sim = 1000
κ_mpc = 1.0e-4

# ## Slow Recovery
obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [Diagonal([12*(t/H_mpc)^2; 2.0*(t/H_mpc)^4]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
	u = [Diagonal([100; 1]) for t = 1:H_mpc-0],
	γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc-0],
	b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

# ## Fast Recovery
obj = TrackingVelocityObjective(model, env, H_mpc,
    q = [Diagonal([12*(t/H_mpc)^2; 12*(t/H_mpc)^2]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
    u = [Diagonal([100; 1]) for t = 1:H_mpc-0],
    γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc-0],
    b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

# ## Policy
p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 10,
		max_time = ref_traj.h/2, # HARD REAL TIME
		),
    mpc_opts = CIMPCOptions());

# ## Disturbances
idx_d1 = 20
idx_d2 = idx_d1 + 200
idx_d3 = idx_d2 + 80
idx_d4 = idx_d3 + 200
idx_d5 = idx_d4 + 30
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
impulses = [[-5.5; 0.0], [+5.5; 0.0], [+5.5; 0.0], [-1.5; 0.0], [-6.5; 0.0]]
d = impulse_disturbances(impulses, idx);

# ## Initial Conditions
q1_sim = [0.0, 0.0]
v1_sim = [0.0; 0.0]

# ## Simulator
sim = simulator(s, H_sim, h=h_sim, policy=p, dist=d)

# ## Simulate
status = simulate!(sim, verbose = true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_robot!(vis, model, sim.traj, sample = 1)
pθ_right = generate_pusher_traj(d, sim.traj, side=:right)
pθ_left  = generate_pusher_traj(d, sim.traj, side=:left)
visualize_disturbance!(vis, model, pθ_right, anim=anim, sample=1, offset=0.05, name=:PusherRight);
visualize_disturbance!(vis, model, pθ_left,  anim=anim, sample=1, offset=0.05, name=:PusherLeft);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio