const ContactControl = Main
vis = Visualizer()
open(vis)

include(joinpath(module_dir(), "src", "dynamics", "walledcartpole", "model.jl"))
include(joinpath(module_dir(), "src", "dynamics", "walledcartpole", "visuals.jl"))

s = get_simulation("walledcartpole", "flat_2D_lc", "flat")
model = s.model
env = s.env

# time
h = 0.04
H = 50

# reference trajectory
ref_traj = contact_trajectory(model, env, H, h)
ref_traj.h
qref = [0.0; 0.0; 0.0; 0.0]
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
	r = residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)#[model.dim.q .+ (1:model.dim.c)]
	@test norm(r) < 1.0e-4
end

# MPC
N_sample = 2
H_mpc = 10
h_sim = h / N_sample
H_sim = 1100

# barrier parameter
κ_mpc = 2.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [[Diagonal([1e-1;1e-3;1e-8;1e-8]) for t = 1:H_mpc-1];
		 [Diagonal([1e+1;1e+0;1e-8;1e-8]) for t = 1:1]],
	v = [[Diagonal([1e-0;3e+1;1e-8;1e-8]) for t = 1:H_mpc-1];
		 [Diagonal([1e+1;1e+1;1e-8;1e-8]) for t = 1:1]],
	u = [Diagonal([3e-2]) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    mpc_opts = LinearizedMPCOptions(),
	n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		max_time = h/5, # we could potentially run the MPC loop 5x faster (5 * 25 = 125 Hz)
		))

# p = open_loop_policy(fill(zeros(nu), H_sim); N_sample = N_sample)

idx_d1 = 20
idx_d2 = idx_d1 + 200
idx_d3 = idx_d2 + 150
idx_d4 = idx_d3 + 200
idx_d5 = idx_d4 + 150
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
impulses = [[+0.2; 0;0;0], [+0.2; 0;0;0], [-0.2; 0;0;0], [+0.2; 0;0;0], [+0.2; 0;0;0]]
d = impulse_disturbances(impulses, idx)

q1_sim = SVector{model.dim.q}([0.00, 0.00, 0.00, 0.00])
q0_sim = SVector{model.dim.q}([0.00, 0.00, 0.00, 0.00])

sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	d = d,
	)

status = ContactControl.simulate!(sim, verbose = true)
anim = visualize_robot!(vis, model, sim.traj, sample = 1)

################################################################################
# Timing result
################################################################################
process!(sim)
# Time budget
ref_traj.h
# Time used on average
sim.stats.μ_dt
# Speed ratio
H_sim * h_sim / sum(sim.stats.dt)



# filename = "walledcartpole_robust"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/video/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/video/$filename.mp4",
#     "/home/simon/Documents/video/$filename.gif", overwrite=true)
