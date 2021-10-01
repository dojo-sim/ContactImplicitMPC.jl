using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays

# Visualizer
vis = ContactImplicitMPC.Visualizer()
open(vis)

# Simulation
s = get_simulation("pushbot", "flat_2D_lc", "flat")
model = s.model
env = s.env

h = 0.04
H = 100

# Reference Trajectory
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h
qref = [0.0; 0.0]
ur = zeros(model.dim.u)
γr = zeros(model.dim.c)
br = zeros(model.dim.c * friction_dim(env))
ψr = zeros(model.dim.c)
ηr = zeros(model.dim.c * friction_dim(env))
wr = zeros(model.dim.w)

# set reference
for t = 1:H
	ref_traj.z[t] = pack_z(model, env, qref, γr, br, ψr, ηr)
	ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
end

# test reference
for t = 1:H
	r = ContactImplicitMPC.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)#[model.dim.q .+ (1:model.dim.c)]
	@assert norm(r) < 1.0e-4
end

# Initial conditions
q0 = @SVector [0.0 * π, 0.0]
q1 = @SVector [0.0 * π, 0.0]

# Simulator
sim = ContactImplicitMPC.simulator(s, q0, q1, h, H,
	ip_opts = ContactImplicitMPC.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-5),
	sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = false))

# Simulate
status = ContactImplicitMPC.simulate!(sim)

# MPC setup 
N_sample = 2
H_mpc = 40
h_sim = h / N_sample
H_sim = 1000

# barrier parameter
κ_mpc = 1.0e-4

# SLOW RECOVERY
obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [Diagonal([12*(t/H_mpc)^2; 2.0*(t/H_mpc)^4]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
	u = [Diagonal([100; 1]) for t = 1:H_mpc-0],
	γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc-0],
	b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

	# FAST RECOVERY
obj = TrackingVelocityObjective(model, env, H_mpc,
    q = [Diagonal([12*(t/H_mpc)^2; 12*(t/H_mpc)^2]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
    u = [Diagonal([100; 1]) for t = 1:H_mpc-0],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc-0],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 10,
		max_time = ref_traj.h/2, # HARD REAL TIME
		# verbose = true,
		),
    mpc_opts = LinearizedMPCOptions(),
    )

# Disturbances
idx_d1 = 20
idx_d2 = idx_d1 + 200
idx_d3 = idx_d2 + 80
idx_d4 = idx_d3 + 200
idx_d5 = idx_d4 + 30
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
impulses = [[-5.5; 0.0], [+5.5; 0.0], [+5.5; 0.0], [-1.5; 0.0], [-6.5; 0.0]]
d = impulse_disturbances(impulses, idx)

# Initial conditions
q1_sim = SVector{model.dim.q}([0.0, 0.0])
q0_sim = SVector{model.dim.q}([0.0, 0.0])

# Simulator
sim = ContactImplicitMPC.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	d = d,
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true))

# Simulate
status = ContactImplicitMPC.simulate!(sim, verbose = true)

################################################################################
# Timing result
################################################################################
process!(sim)
# Time budget
ref_traj.h # 0.04
# Time used on average
sim.stats.μ_dt # 0.0138
sim.stats.σ_dt # 0.0138
# Speed ratio
H_sim * h_sim / sum(sim.stats.dt) # 2.90

# Visualize
anim = visualize_robot!(vis, model, sim.traj, sample = 1)
pθ_right = generate_pusher_traj(d, sim.traj, side=:right)
pθ_left  = generate_pusher_traj(d, sim.traj, side=:left)
visualize_disturbance!(vis, model, pθ_right, anim=anim, sample=1, offset=0.05, name=:PusherRight)
visualize_disturbance!(vis, model, pθ_left,  anim=anim, sample=1, offset=0.05, name=:PusherLeft)