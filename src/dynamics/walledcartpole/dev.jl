const ContactControl = Main
vis = Visualizer()
open(vis)

include(joinpath(module_dir(), "src", "dynamics", "walledcartpole", "model.jl"))
include(joinpath(module_dir(), "src", "dynamics", "walledcartpole", "visuals.jl"))

s = get_simulation("walledcartpole", "flat_2D_lc", "flat")
model = s.model
env = s.env

# time
h = 0.01
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
J_func(model,env, q0)
# MPC
N_sample = 2
H_mpc = 10
h_sim = h / N_sample
H_sim = 1500

# barrier parameter
κ_mpc = 2.0e-4

# SLOW RECOVERY
obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [[Diagonal([1; 1e-4]) for t = 1:H_mpc-1];
		 [Diagonal([10; 1]) for t = 1:1]],
	v = [[Diagonal([1e-2; 1e-4] ./ (h^2.0)) for t = 1:H_mpc-1];
		 [Diagonal([1e-1; 1e-1] ./ (h^2.0)) for t = 1:1]],
	u = [[Diagonal([5e-3]) for t = 1:H_mpc-1];
		 [Diagonal([5e-2]) for t = 1:1]])

# p = linearized_mpc_policy(ref_traj, s, obj,
#     H_mpc = H_mpc,
#     N_sample = N_sample,
#     κ_mpc = κ_mpc,
#     mpc_opts = LinearizedMPCOptions())

p = open_loop_policy(fill(zeros(nu), H_sim); N_sample = N_sample)

idx_d1 = 20
idx_d2 = idx_d1 + 1500
idx_d3 = idx_d2 + 650
idx_d4 = idx_d3 + 600
idx_d5 = idx_d4 + 650
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
impulses = [[+0.1, 0.0], [+0.1; 0.0], [+0.1; 0.0], [+0.12; 0.0], [+0.12; 0.0]]
d = impulse_disturbances(impulses, idx)

q1_sim = SVector{model.dim.q}([0.00, 0.00, 0.00, 0.00])
q0_sim = SVector{model.dim.q}([0.00, 0.00, 0.00, 0.00])
# q1_sim = SVector{model.dim.q}([0.00, 0.2])
# q0_sim = SVector{model.dim.q}([0.02, 0.18])


sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	d = d
	)

telap = @elapsed status = ContactControl.simulate!(sim, verbose = true)

anim = visualize_robot!(vis, model, sim.traj, sample = 1)


l = 1
lu = 1
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[l:l] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][lu:lu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[3,1], hcat(Vector.(vcat([fill(ref_traj.γ[i][1:nc], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[lu:lu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)

sample = 1

anim = visualize_robot!(vis, model, sim.traj, sample = 1)#, sample = sample, α=1.0)
# pθ_right = generate_pusher_traj(d, sim.traj, side=:right)
# pθ_left  = generate_pusher_traj(d, sim.traj, side=:left)
# visualize_disturbance!(vis, model, pθ_right, anim=anim, sample=sample, offset=0.05, name=:PusherRight)
# visualize_disturbance!(vis, model, pθ_left,  anim=anim, sample=sample, offset=0.05, name=:PusherLeft)



# filename = "walledcartpole_sim"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/video/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/video/$filename.mp4",
#     "/home/simon/Documents/video/$filename.gif", overwrite=true)
